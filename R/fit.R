# Fitting the model to daily hospital admissions.
#
# The odin models compare modelled daily admissions with data through a
# negative binomial likelihood (see the end of inst/odin/*.R). dust2 turns
# that into a likelihood (deterministic "unfilter" or stochastic particle
# filter), and monty provides priors and samplers.

# Parameters that can be fitted. Each entry describes
#   prior:       default prior (monty DSL distribution)
#   description: shown by fit_parameters()
#   initial:     starting value, function(model, nb_size)
#   apply:       backend parameters it sets, function(value, ctx) -> list
#                (list() if the fitted name is itself a backend parameter)
#   rebuild:     build_model() arguments at a given value, function(args, value)
# To add a fittable parameter, add an entry here (see
# docs/modifying_the_model.md).
fit_param_specs <- list(
  R0 = list(
    prior = "Uniform(0.5, 6)",
    description = "Basic reproduction number (beta is recalibrated from it)",
    initial = function(model, nb_size) R0(model),
    apply = function(value, ctx) list(beta = value * ctx$beta_per_R0),
    rebuild = function(args, value) {
      args$R0 <- value
      args$beta <- NULL
      args
    }),
  hosp_scale = list(
    prior = "Uniform(0.1, 10)",
    description = "Multiplier on p_hosp for all groups (capped at 1)",
    initial = function(model, nb_size) 1,
    apply = function(value, ctx) list(p_hosp = pmin(ctx$p_hosp * value, 1)),
    rebuild = function(args, value) {
      args$p_hosp <- pmin(args$p_hosp * value, 1)
      args
    }),
  nb_size = list(
    prior = "Uniform(1, 200)",
    description = "Negative binomial size (dispersion) of the admissions data",
    initial = function(model, nb_size) nb_size,
    # same name as the backend parameter, so the packer passes it on directly
    apply = function(value, ctx) list(),
    rebuild = function(args, value) args))

#' Parameters that can be fitted
#'
#' @return A data.frame with the parameter names that can be passed to
#'   `fit_model(pars = ...)`, their default priors and a description.
#' @export
fit_parameters <- function() {
  data.frame(
    name = names(fit_param_specs),
    default_prior = vapply(fit_param_specs, `[[`, "", "prior"),
    description = vapply(fit_param_specs, `[[`, "", "description"),
    row.names = NULL)
}

#' Fit a model to daily hospital admissions
#'
#' Estimates `R0` (and optionally other parameters, see [fit_parameters()])
#' from observed daily hospital admissions, keeping all other parameters of
#' the model fixed. The observations are compared with modelled admissions
#' through a negative binomial distribution.
#'
#' Three methods are available:
#'
#' * `"optim"`: maximum a posteriori point estimate from the deterministic
#'   model, with approximate 95% intervals from the curvature of the
#'   log posterior. Fast; good for a first look.
#' * `"mcmc"`: Bayesian MCMC on the deterministic model with monty's adaptive
#'   random-walk sampler, started from the `"optim"` estimate. Gives
#'   posterior distributions.
#' * `"pmcmc"`: particle MCMC on the stochastic model: each likelihood is
#'   estimated with a particle filter, so it also accounts for randomness in
#'   the epidemic itself. Much slower; the proposal is tuned from the
#'   deterministic fit.
#'
#' @param model A model from [build_model()]. Its parameters are the fixed
#'   values for everything that is not fitted, and its time 0 is day 0 of the
#'   data.
#' @param data A data.frame with columns `time` (whole days since the start
#'   of the model, >= 1) and `admissions` (new admissions during day
#'   `(time - 1, time]`; `NA` = missing). Optional `region` and/or `age`
#'   columns split the data into separate series, e.g. admissions by region;
#'   without them the data are national totals.
#' @param pars Parameters to fit, see [fit_parameters()].
#' @param method `"optim"`, `"mcmc"` or `"pmcmc"`, see Details.
#' @param prior Optional priors replacing the defaults: a named character
#'   vector of monty distributions, e.g. `c(R0 = "Normal(2, 0.5)")`, or a
#'   complete prior from [monty::monty_dsl()] over `pars`.
#' @param initial Optional named starting values (default: the model's value
#'   for `R0`, 1 for `hosp_scale`, `nb_size` for `nb_size`).
#' @param nb_size Negative binomial size used when `nb_size` is not fitted.
#' @param n_steps,n_chains,burnin Number of MCMC steps per chain, number of
#'   chains, and number of initial steps discarded from each chain (default:
#'   a quarter of `n_steps`).
#' @param n_particles,dt,n_threads Particle filter settings for `"pmcmc"`:
#'   number of particles, time step of the stochastic model, and threads.
#' @param seed Random seed (sets R's seed; `NULL` = leave it alone).
#'
#' @return An object of class `seirvax_fit`. `fit$model` is the model rebuilt
#'   at the estimates (posterior median for MCMC), ready for [run_model()].
#'   See also `summary()`, `plot()` and [project_fit()].
#' @export
fit_model <- function(model, data, pars = "R0",
                      method = c("optim", "mcmc", "pmcmc"),
                      prior = NULL, initial = NULL, nb_size = 10,
                      n_steps = 2000, n_chains = 4, burnin = NULL,
                      n_particles = 200, dt = 0.1, n_threads = 1,
                      seed = NULL) {
  if (!inherits(model, "seirvax_model")) {
    stop("'model' must be created with build_model()")
  }
  method <- match.arg(method)
  pars <- unique(pars)
  unknown <- setdiff(pars, names(fit_param_specs))
  if (length(pars) == 0 || length(unknown) > 0) {
    stop(sprintf("'pars' must be one or more of: %s",
                 paste(names(fit_param_specs), collapse = ", ")))
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  obs <- prepare_fit_data(data, model)
  posterior_det <- fit_posterior(model, obs, pars, prior, nb_size,
                                 type = "deterministic")
  x0 <- fit_initial(model, pars, initial, nb_size, posterior_det$posterior)

  domain <- posterior_det$posterior$domain
  opt <- fit_optim(posterior_det$posterior, x0)
  ret <- list(method = method, pars = pars, data = obs$data_long,
              series = obs$series, base_model = model, nb_size = nb_size,
              domain = domain, optim = opt, samples = NULL)

  if (method == "optim") {
    est <- opt$par
    ret$estimates <- data.frame(
      parameter = pars, estimate = opt$par,
      lower = pmax(opt$par - 1.96 * sqrt(diag(opt$vcv)), domain[, 1]),
      upper = pmin(opt$par + 1.96 * sqrt(diag(opt$vcv)), domain[, 2]),
      row.names = NULL)
  } else {
    burnin <- if (is.null(burnin)) floor(n_steps / 4) else burnin
    if (burnin >= n_steps) {
      stop("'burnin' must be smaller than 'n_steps'")
    }
    if (method == "mcmc") {
      post <- posterior_det$posterior
      sampler <- monty::monty_sampler_adaptive(initial_vcv = opt$vcv)
    } else {
      post <- fit_posterior(model, obs, pars, prior, nb_size,
                            type = "stochastic", n_particles = n_particles,
                            dt = dt, n_threads = n_threads)$posterior
      sampler <- monty::monty_sampler_random_walk(vcv = opt$vcv)
    }
    samples <- monty::monty_sample(post, sampler, n_steps,
                                   initial = opt$par, n_chains = n_chains,
                                   burnin = burnin)
    x <- matrix(samples$pars, nrow = length(pars))
    est <- apply(x, 1, stats::median)
    ret$samples <- samples
    ret$estimates <- data.frame(
      parameter = pars, estimate = est,
      lower = apply(x, 1, stats::quantile, 0.025, names = FALSE),
      upper = apply(x, 1, stats::quantile, 0.975, names = FALSE),
      row.names = NULL)
    density <- matrix(samples$density, ncol = n_chains)
    ret$acceptance <- mean(apply(density, 2, function(d) mean(diff(d) != 0)))
    ret$n_particles <- if (method == "pmcmc") n_particles else NULL
    ret$dt <- if (method == "pmcmc") dt else NULL
  }
  names(est) <- pars
  ret$model <- fit_rebuild(model, est)
  if ("nb_size" %in% pars) {
    ret$nb_size <- est[["nb_size"]]
  }
  structure(ret, class = "seirvax_fit")
}

# Turn the long data.frame into dust2's format: one row per time, with the
# admissions of the n_obs series as a list column; plus the obs_map matrix.
prepare_fit_data <- function(data, model) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!all(c("time", "admissions") %in% names(data))) {
    stop("'data' must have columns 'time' and 'admissions'")
  }
  extra <- setdiff(names(data), c("time", "admissions", "age", "region"))
  if (length(extra) > 0) {
    stop(sprintf("Unknown column(s) in 'data': %s", paste(extra, collapse = ", ")))
  }
  if (anyNA(data$time) || any(data$time < 1) ||
      any(data$time != round(data$time))) {
    stop("'data$time' must be whole days >= 1 (day 0 is the model start)")
  }
  a <- data$admissions[!is.na(data$admissions)]
  if (any(a < 0) || any(a != round(a))) {
    stop("'data$admissions' must be non-negative whole numbers (or NA)")
  }

  dims <- model$dims
  keys <- intersect(c("age", "region"), names(data))
  for (k in keys) {
    data[[k]] <- as.character(data[[k]])
    bad <- setdiff(stats::na.omit(data[[k]]), dims[[k]])
    if (length(bad) > 0) {
      stop(sprintf("Unknown %s in 'data': %s", k, paste(bad, collapse = ", ")))
    }
  }
  if (length(keys) == 0) {
    series <- data.frame(label = "all")
    data$.series <- 1L
  } else {
    series <- unique(data[keys])
    rownames(series) <- NULL
    key <- do.call(paste, c(data[keys], sep = "\r"))
    data$.series <- match(key, do.call(paste, c(series, sep = "\r")))
    series$label <- apply(series[keys], 1, function(x) {
      paste(ifelse(is.na(x), "all", x), collapse = " / ")
    })
  }
  if (anyDuplicated(data[c("time", ".series")])) {
    stop("'data' has more than one row for the same time and series")
  }

  n_obs <- nrow(series)
  group_age <- rep(dims$age, length(dims$region))
  group_region <- rep(dims$region, each = length(dims$age))
  obs_map <- matrix(0, n_obs, length(group_age))
  for (k in seq_len(n_obs)) {
    match <- rep(TRUE, length(group_age))
    if (!is.null(series$age) && !is.na(series$age[k])) {
      match <- match & group_age == series$age[k]
    }
    if (!is.null(series$region) && !is.na(series$region[k])) {
      match <- match & group_region == series$region[k]
    }
    obs_map[k, match] <- 1
  }

  times <- sort(unique(data$time))
  y <- matrix(NA_real_, length(times), n_obs)
  y[cbind(match(data$time, times), data$.series)] <- data$admissions
  data_dust <- data.frame(time = times)
  data_dust$admissions <- I(lapply(seq_along(times), function(i) y[i, ]))

  data$series <- series$label[data$.series]
  data$.series <- NULL
  list(data_dust = data_dust, obs_map = obs_map, series = series,
       data_long = data)
}

fit_posterior <- function(model, obs, pars, prior, nb_size, type,
                          n_particles = NULL, dt = NULL, n_threads = 1) {
  fixed <- model$pars
  fixed$n_obs <- nrow(obs$obs_map)
  fixed$obs_map <- obs$obs_map
  fixed$nb_size <- nb_size
  ctx <- list(beta_per_R0 = model$beta / R0(model), p_hosp = fixed$p_hosp)
  specs <- fit_param_specs[pars]

  # Backend parameters set by the fitted ones must not also be fixed
  probe <- unlist(lapply(specs, function(s) names(s$apply(1, ctx))))
  fixed[c(pars, probe)] <- NULL
  packer <- monty::monty_packer(pars, fixed = fixed, process = function(x) {
    unlist(lapply(pars, function(p) specs[[p]]$apply(x[[p]], ctx)),
           recursive = FALSE)
  })

  if (type == "deterministic") {
    filter <- dust2::dust_unfilter_create(seir_ode(), 0, obs$data_dust)
  } else {
    filter <- dust2::dust_filter_create(seir_stoch(), 0, obs$data_dust,
                                        n_particles = n_particles, dt = dt,
                                        n_threads = n_threads)
  }
  likelihood <- dust2::dust_likelihood_monty(filter, packer)
  list(posterior = likelihood + fit_prior(pars, prior), packer = packer)
}

fit_prior <- function(pars, prior) {
  if (inherits(prior, "monty_model")) {
    if (!identical(prior$parameters, pars)) {
      stop(sprintf("The prior must be over exactly these parameters, in order: %s",
                   paste(pars, collapse = ", ")))
    }
    return(prior)
  }
  dist <- vapply(fit_param_specs[pars], `[[`, "", "prior")
  if (!is.null(prior)) {
    prior <- unlist(prior)
    unknown <- setdiff(names(prior), pars)
    if (is.null(names(prior)) || length(unknown) > 0) {
      stop("'prior' must be named with parameters being fitted")
    }
    dist[names(prior)] <- prior
  }
  # monty_dsl() captures its argument unevaluated; do.call passes the string
  do.call(monty::monty_dsl, list(paste(pars, "~", dist, collapse = "\n")))
}

fit_initial <- function(model, pars, initial, nb_size, posterior) {
  x0 <- vapply(pars, function(p) fit_param_specs[[p]]$initial(model, nb_size), 1)
  if (!is.null(initial)) {
    unknown <- setdiff(names(initial), pars)
    if (is.null(names(initial)) || length(unknown) > 0) {
      stop("'initial' must be named with parameters being fitted")
    }
    x0[names(initial)] <- unlist(initial)
  }
  if (!is.finite(monty::monty_model_density(posterior, x0))) {
    stop("The starting values have zero prior probability or likelihood; ",
         "set 'initial' or widen the prior")
  }
  x0
}

# Maximum a posteriori estimate and a proposal / uncertainty VCV from the
# curvature at the optimum.
fit_optim <- function(posterior, x0) {
  f <- function(x) {
    v <- monty::monty_model_density(posterior, x)
    if (is.finite(v)) v else -1e10
  }
  dom <- posterior$domain
  opt <- stats::optim(x0, f, method = "L-BFGS-B",
                      lower = dom[, 1], upper = dom[, 2],
                      control = list(fnscale = -1, parscale = pmax(abs(x0), 0.1)),
                      hessian = TRUE)
  vcv <- tryCatch(solve(-opt$hessian), error = function(e) NULL)
  if (is.null(vcv) || any(!is.finite(vcv)) ||
      any(eigen(vcv, only.values = TRUE)$values <= 0)) {
    warning("Could not estimate uncertainty from the curvature at the optimum ",
            "(is a parameter at the edge of its prior?); using a rough guess")
    vcv <- diag((0.05 * pmax(abs(opt$par), 0.1))^2, length(x0))
  }
  names(opt$par) <- names(x0)
  dimnames(vcv) <- list(names(x0), names(x0))
  list(par = opt$par, log_posterior = opt$value, vcv = vcv,
       convergence = opt$convergence, message = opt$message)
}

fit_rebuild <- function(model, values) {
  args <- model$args
  for (p in names(values)) {
    args <- fit_param_specs[[p]]$rebuild(args, values[[p]])
  }
  do.call(build_model, args)
}

# Parameter draws: posterior samples for MCMC; for optim, draws from the
# normal approximation at the estimate, restricted to the prior's domain.
fit_draws <- function(fit, n_draws) {
  if (is.null(fit$samples)) {
    L <- chol(fit$optim$vcv)
    x <- matrix(0, 0, length(fit$pars))
    for (i in 1:100) {
      z <- matrix(stats::rnorm(n_draws * length(fit$pars)), n_draws)
      y <- sweep(z %*% L, 2, fit$optim$par, "+")
      ok <- apply(y, 1, function(v) all(v > fit$domain[, 1] & v < fit$domain[, 2]))
      x <- rbind(x, y[ok, , drop = FALSE])
      if (nrow(x) >= n_draws) break
    }
    x <- x[seq_len(min(n_draws, nrow(x))), , drop = FALSE]
  } else {
    x <- t(matrix(fit$samples$pars, nrow = length(fit$pars)))
    x <- x[sample.int(nrow(x), min(n_draws, nrow(x))), , drop = FALSE]
  }
  colnames(x) <- fit$pars
  x
}

#' Project a fitted model
#'
#' Runs the model for parameter draws from the fit (posterior samples for
#' MCMC fits; draws from the normal approximation at the estimate for
#' `method = "optim"`), e.g. to project beyond the data. Each
#' draw restarts from time 0; for `"pmcmc"` fits the projections are
#' therefore not conditioned on the particular trajectory the data followed.
#'
#' @param fit Result of [fit_model()].
#' @param t_max Last day to run to (default: last day of the data).
#' @param n_draws Number of posterior draws (MCMC fits).
#' @param type `"deterministic"` or `"stochastic"`.
#' @param n_sim Stochastic simulations per draw.
#' @param vars,by Output to keep, as in [get_output()].
#' @param seed Random seed for choosing the draws (`NULL` = leave R's seed
#'   alone, so repeated calls give slightly different bands).
#' @param ... Passed to [run_model()] (e.g. `dt`, `n_threads`).
#'
#' @return A `data.table` as from [get_output()] with an extra `draw` column.
#'   The parameter values behind each `draw` are attached as
#'   `attr(x, "draws")`.
#' @export
project_fit <- function(fit, t_max = NULL, n_draws = 100,
                        type = c("deterministic", "stochastic"), n_sim = 1,
                        vars = "inc_hosp", by = character(0), seed = NULL,
                        ...) {
  type <- match.arg(type)
  if (is.null(t_max)) {
    t_max <- max(fit$data$time)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  draws <- fit_draws(fit, n_draws)
  res <- lapply(seq_len(nrow(draws)), function(i) {
    m <- fit_rebuild(fit$base_model, draws[i, ])
    r <- run_model(m, type, t_max = t_max, n_sim = n_sim, vars = vars, ...)
    get_output(r, vars, by)[, draw := i]
  })
  res <- rbindlist(res)
  setcolorder(res, c("time", "draw"))
  # keep the parameter values, so callers can pair a trajectory with the
  # parameters that produced it (plot() needs nb_size per draw)
  setattr(res, "draws", draws)
  res[]
}

#' Outcomes of a fitted model, by age and region
#'
#' Totals at the end of the projection (infections, hospital admissions and
#' deaths) and the peak of hospital occupancy, broken down by age group
#' and/or region, with the uncertainty from the fit. Each posterior draw gives
#' one set of outcomes; the table reports the median and a central interval
#' over draws.
#'
#' The `share_*` columns give each stratum's share of the national total,
#' computed within each draw and then summarised, so they answer "what
#' fraction of admissions are in 65+, and how sure are we?" rather than
#' dividing two uncertain totals. Shares are reported for the totals only, not
#' for the peak measures: strata peak at different times, so their peaks do
#' not add up to the national peak.
#'
#' @inheritParams project_fit
#' @param by Stratifications to report, any of `"age"`, `"region"` and
#'   `"vacc_status"`. The default, `character(0)`, gives national totals.
#' @param measures Which outcomes to report (default: all). One or more of
#'   `"cum_infections"`, `"cum_hosp"` (admissions), `"D"` (deaths),
#'   `"peak_hosp_occupancy"` and `"peak_hosp_day"`. They are named after the
#'   output variables of [get_output()].
#' @param level Width of the interval.
#'
#' @return A `data.table` with the `by` columns, `measure`, `mean`, `median`,
#'   `lower` and `upper`, and the same three summaries of `share` when `by` is
#'   given.
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_model(model, data, pars = "R0", method = "mcmc")
#' project_summary(fit, by = "age", t_max = 365)
#' project_summary(fit, by = c("region", "age"), measures = "cum_hosp")
#' }
project_summary <- function(fit, by = character(0), t_max = NULL,
                            n_draws = 100,
                            type = c("deterministic", "stochastic"),
                            n_sim = 1, measures = NULL, level = 0.95,
                            seed = NULL, ...) {
  by <- check_by(by)
  if (is.null(measures)) {
    measures <- outcome_measures
  }
  unknown <- setdiff(measures, outcome_measures)
  if (length(unknown) > 0) {
    stop(sprintf("'measures' must be a subset of: %s",
                 paste(outcome_measures, collapse = ", ")))
  }
  vars <- c("cum_infections", "cum_hosp", "D", "H")
  proj <- project_fit(fit, t_max = t_max, n_draws = n_draws, type = type,
                      n_sim = n_sim, vars = vars, by = by, seed = seed, ...)
  per_draw <- trajectory_measures(proj, c("draw", "sim"), by)
  per_draw <- per_draw[measure %in% measures]

  ret <- summarise_measures(per_draw, by, level)
  # Shares only make sense for the totals: strata peak at different times, so
  # a "share of the peak" would not add up to anything meaningful.
  share_of <- intersect(measures, c("cum_infections", "cum_hosp", "D"))
  if (length(by) > 0 && length(share_of) > 0) {
    sh <- per_draw[measure %in% share_of]
    sh[, value := value / sum(value), by = c("draw", "sim", "measure")]
    shares <- summarise_measures(sh, by, level)
    setnames(shares, c("mean", "median", "lower", "upper"),
             paste0("share_", c("mean", "median", "lower", "upper")))
    ret <- merge(ret, shares, by = c(by, "measure"), all.x = TRUE)
  }
  setorderv(ret, c("measure", by))
  ret[, measure := as.character(measure)]
  ret[]
}

#' @export
print.seirvax_fit <- function(x, ...) {
  cat(sprintf("<seirvax_fit> method: %s\n", x$method))
  cat(sprintf("  Data: %d day(s) of admissions, %d series (%s)\n",
              length(unique(x$data$time)), nrow(x$series),
              paste(x$series$label, collapse = ", ")))
  if (!is.null(x$samples)) {
    d <- dim(x$samples$pars)
    cat(sprintf("  MCMC: %d chain(s) x %d step(s) after burn-in; acceptance %.2f\n",
                d[3], d[2], x$acceptance))
    if (x$method == "pmcmc") {
      cat(sprintf("  Particle filter: %d particles, dt = %g\n",
                  x$n_particles, x$dt))
    }
  } else {
    cat(sprintf("  Log posterior at estimate: %.2f\n", x$optim$log_posterior))
  }
  est <- x$estimates
  est[-1] <- lapply(est[-1], signif, 4)
  cat(sprintf("  Estimates (%s):\n", if (is.null(x$samples))
    "maximum a posteriori, approx. 95% interval" else
    "posterior median, 95% credible interval"))
  print(est, row.names = FALSE)
  cat("  Fitted model: fit$model; see also plot(fit) and project_fit()\n")
  invisible(x)
}

#' @export
summary.seirvax_fit <- function(object, ...) {
  object$estimates
}

#' Plot a fit
#'
#' `type = "fit"` shows the data with the modelled admissions: the median and
#' 95% interval of the expected admissions over parameter draws (see
#' [project_fit()]), and a lighter 95% predictive interval that also
#' includes the negative binomial observation noise. `type = "trace"` shows
#' MCMC traces. `type = "distribution"` shows one outcome broken down by age
#' and/or region, with its uncertainty (see [project_summary()]).
#'
#' @param x Result of [fit_model()].
#' @param type `"fit"`, `"trace"` or `"distribution"`.
#' @param by Stratification for `type = "distribution"`: the first is the
#'   axis, a second is shown by colour.
#' @param measure Outcome for `type = "distribution"`, see
#'   [project_summary()].
#' @param t_max Last day to show (default: last day of the data).
#' @param n_draws Number of posterior draws used.
#' @param seed Random seed for the draws and the observation noise.
#' @param ... Unused.
#' @return A ggplot object.
#' @export
plot.seirvax_fit <- function(x, type = c("fit", "trace", "distribution"),
                             t_max = NULL, n_draws = 100, by = "age",
                             measure = "cum_hosp", seed = NULL, ...) {
  type <- match.arg(type)
  if (type == "distribution") {
    if (length(measure) != 1) {
      stop("'measure' must be a single outcome")
    }
    by <- check_by(by)
    if (length(by) == 0) {
      stop("'by' is needed for a distribution plot, e.g. by = \"age\"")
    }
    d <- project_summary(x, by = by, t_max = t_max, n_draws = n_draws,
                         measures = measure, seed = seed)
    # Every aesthetic belongs to the pointrange layer, and none to the plot,
    # so observed data can be added as another layer without inheriting
    # mappings its own data cannot supply.
    mapping <- ggplot2::aes(x = .data[[by[[1]]]], y = median,
                            ymin = lower, ymax = upper)
    p <- ggplot2::ggplot(d)
    if (length(by) > 1) {
      mapping <- utils::modifyList(
        mapping, ggplot2::aes(colour = .data[[by[[2]]]]))
      p <- p + ggplot2::geom_pointrange(
        mapping, position = ggplot2::position_dodge(width = 0.4))
    } else {
      p <- p + ggplot2::geom_pointrange(mapping)
    }
    return(p + ggplot2::labs(x = NULL, y = measure) +
             ggplot2::expand_limits(y = 0) + ggplot2::theme_bw())
  }
  if (type == "trace") {
    if (is.null(x$samples)) {
      stop("No MCMC samples: fit with method = \"mcmc\" or \"pmcmc\"")
    }
    p <- x$samples$pars
    d <- data.table(parameter = rep(x$pars, times = prod(dim(p)[2:3])),
                    step = rep(seq_len(dim(p)[2]), each = dim(p)[1]),
                    chain = factor(rep(seq_len(dim(p)[3]), each = prod(dim(p)[1:2]))),
                    value = as.vector(p))
    return(ggplot2::ggplot(d, ggplot2::aes(step, value, colour = chain)) +
             ggplot2::geom_line(alpha = 0.7) +
             ggplot2::facet_wrap(~parameter, scales = "free_y") +
             ggplot2::labs(x = "Step (after burn-in)", y = NULL) +
             ggplot2::theme_bw())
  }

  proj <- project_fit(x, t_max = t_max, n_draws = n_draws,
                      by = c("age", "region"), seed = seed)
  draws <- attr(proj, "draws")   # the draws behind the trajectories
  nb_size <- if ("nb_size" %in% x$pars) NULL else x$nb_size
  series <- x$series
  mod <- rbindlist(lapply(seq_len(nrow(series)), function(k) {
    sel <- rep(TRUE, nrow(proj))
    if (!is.null(series$age) && !is.na(series$age[k])) {
      sel <- sel & proj$age == series$age[k]
    }
    if (!is.null(series$region) && !is.na(series$region[k])) {
      sel <- sel & proj$region == series$region[k]
    }
    proj[sel, .(value = sum(value)), by = .(time, draw)][, series := series$label[k]]
  }))
  mod <- mod[time > 0]
  # several observation-noise samples per draw for a smooth predictive band
  mod <- mod[rep(seq_len(.N), ceiling(1000 / nrow(draws)))]
  size <- if (is.null(nb_size)) draws[mod$draw, "nb_size"] else nb_size
  mod[, predicted := stats::rnbinom(.N, size = size, mu = value)]
  q <- mod[, .(median = stats::median(value),
               lower = stats::quantile(value, 0.025, names = FALSE),
               upper = stats::quantile(value, 0.975, names = FALSE),
               pred_lower = stats::quantile(predicted, 0.025, names = FALSE),
               pred_upper = stats::quantile(predicted, 0.975, names = FALSE)),
           by = .(series, time)]
  obs <- as.data.table(x$data)[!is.na(admissions)]
  series_levels <- x$series$label
  q[, series := factor(series, series_levels)]
  obs[, series := factor(series, series_levels)]

  p <- ggplot2::ggplot(q, ggplot2::aes(time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_lower, ymax = pred_upper),
                         fill = "steelblue", alpha = 0.15) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "steelblue", alpha = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = median), colour = "steelblue4") +
    ggplot2::geom_point(ggplot2::aes(y = admissions), data = obs, size = 1) +
    ggplot2::labs(x = "Day", y = "Daily hospital admissions") +
    ggplot2::theme_bw()
  if (nrow(series) > 1) {
    p <- p + ggplot2::facet_wrap(~series, scales = "free_y")
  }
  p
}
