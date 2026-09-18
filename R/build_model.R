#' Build a model
#'
#' Collects, checks and expands all parameters of the SEIR-HD model with a
#' vaccinated stratum. Nothing is simulated; pass the result to [run_model()].
#'
#' The population is stratified by age group and region. Most parameters can
#' be given as a single value, a vector by age group, a vector by region
#' (length `n_region`, or named with the region labels) or an
#' `n_age x n_region` matrix.
#'
#' Vaccination is a fixed initial coverage. Vaccinated people differ from
#' unvaccinated people only through the vaccine effectiveness (VE) arguments,
#' which can also vary by age and region:
#'
#' * `ve_infection`: susceptibility is multiplied by `1 - ve_infection`.
#' * `ve_transmission`: infectiousness of vaccinated infected people is
#'   multiplied by `1 - ve_transmission`.
#' * `ve_hosp`: `p_hosp` is multiplied by `1 - ve_hosp`. This is protection
#'   *given infection*; the total protection against hospitalisation is
#'   `1 - (1 - ve_infection) * (1 - ve_hosp)`.
#' * `ve_death`: `p_death_hosp` and `p_death_comm` are multiplied by
#'   `1 - ve_death` (again given infection / hospitalisation).
#'
#' @param N Population: a single number, a (named) vector by age group, or an
#'   `n_age x n_region` matrix whose dimnames give the age group and region
#'   labels. A region-only model is a `1 x n_region` matrix.
#' @param C_age Age contact matrix, `C_age[a, b]` = contacts per day an
#'   age-`a` person has with age-`b` people. Required with more than one age
#'   group. Only its shape matters when `R0` is given, not its scale.
#' @param M_region Region mixing matrix, `M_region[r, s]` = fraction of a
#'   region-`r` resident's contacts made with region-`s` residents (rows sum to
#'   1). Required with more than one region. See [combine_mixing()].
#' @param R0,beta Give exactly one. `R0` calibrates `beta` with the
#'   next-generation matrix (see [R0()]).
#' @param latent_period,infectious_period Mean durations (days) of E and I.
#' @param susceptibility Relative susceptibility of unvaccinated people.
#' @param transmissibility Relative infectiousness of unvaccinated infected
#'   people.
#' @param p_hosp Probability that an unvaccinated infection is hospitalised.
#' @param p_death_hosp Probability that an unvaccinated hospitalised case dies.
#' @param p_death_comm Probability that an unvaccinated infection that is not
#'   hospitalised dies (in the community).
#' @param dur_hosp Mean length of hospital stay (days).
#' @param dur_hosp_vacc Mean length of stay for vaccinated patients (default:
#'   `dur_hosp`).
#' @param vacc_coverage Fraction of each group vaccinated at the start.
#' @param ve_infection,ve_transmission,ve_hosp,ve_death Vaccine effectiveness,
#'   see Details.
#' @param I0,E0 Initial number of infectious / exposed people (placed in the
#'   unvaccinated stratum).
#' @param beta_schedule Optional contact multipliers over time, e.g. for
#'   interventions. Either a data.frame with columns `time` and `mult` (applies
#'   to everyone), or with extra columns `age` and/or `region` to target
#'   groups (`NA` = all). Rows are applied in time order: from `time` onwards
#'   the matching groups get multiplier `mult` until a later row matching them
#'   overrides it. `mult` is the fraction of the group's contacts that still
#'   happen: contacts between groups i and j are scaled by
#'   `min(mult_i, mult_j)`, so a group with multiplier m has all its contacts
#'   (with any group) reduced by at least that much.
#'
#' @return An object of class `seirvax_model`.
#' @export
#' @examples
#' m <- build_model(N = 1e5, R0 = 2, I0 = 10, p_hosp = 0.02)
#' m
build_model <- function(N,
                        C_age = NULL,
                        M_region = NULL,
                        R0 = NULL,
                        beta = NULL,
                        latent_period = 3,
                        infectious_period = 5,
                        susceptibility = 1,
                        transmissibility = 1,
                        p_hosp = 0,
                        p_death_hosp = 0,
                        p_death_comm = 0,
                        dur_hosp = 7,
                        dur_hosp_vacc = NULL,
                        vacc_coverage = 0,
                        ve_infection = 0,
                        ve_transmission = 0,
                        ve_hosp = 0,
                        ve_death = 0,
                        I0 = 0,
                        E0 = 0,
                        beta_schedule = NULL) {
  args <- as.list(environment())

  pop <- parse_population(N)
  N <- pop$N
  dims <- pop$dims
  n_group <- length(N)

  C_age <- check_C_age(C_age, length(dims$age), dims$age)
  M_region <- check_M_region(M_region, length(dims$region), dims$region)
  contact <- combine_mixing(C_age, M_region)

  if (is.null(R0) == is.null(beta)) {
    stop("Give exactly one of 'R0' and 'beta'")
  }
  check_range(latent_period, "latent_period", min_open = TRUE)
  check_range(infectious_period, "infectious_period", min_open = TRUE)

  # Everything group-specific as an age x region matrix
  g <- function(x, name, ...) {
    ret <- expand_param(x, name, dims)
    check_range(ret, name, ...)
    ret
  }
  susceptibility <- g(susceptibility, "susceptibility")
  transmissibility <- g(transmissibility, "transmissibility")
  p_hosp <- g(p_hosp, "p_hosp", max = 1)
  p_death_hosp <- g(p_death_hosp, "p_death_hosp", max = 1)
  p_death_comm <- g(p_death_comm, "p_death_comm", max = 1)
  dur_hosp <- g(dur_hosp, "dur_hosp", min_open = TRUE)
  dur_hosp_vacc <- if (is.null(dur_hosp_vacc)) dur_hosp else
    g(dur_hosp_vacc, "dur_hosp_vacc", min_open = TRUE)
  vacc_coverage <- g(vacc_coverage, "vacc_coverage", max = 1)
  ve_infection <- g(ve_infection, "ve_infection", max = 1)
  ve_transmission <- g(ve_transmission, "ve_transmission", max = 1)
  ve_hosp <- g(ve_hosp, "ve_hosp", max = 1)
  ve_death <- g(ve_death, "ve_death", max = 1)
  I0 <- g(I0, "I0")
  E0 <- g(E0, "E0")
  if (any(abs(c(I0, E0) - round(c(I0, E0))) > 1e-8)) {
    stop("'I0' and 'E0' must be whole numbers of people")
  }

  # [group, vaccination status] matrices for the backend
  vs <- function(unvacc, vacc) cbind(as.vector(unvacc), as.vector(vacc))

  S_vacc <- round(N * vacc_coverage)
  S_unvacc <- N - S_vacc - I0 - E0
  if (any(S_unvacc < 0)) {
    stop("I0 + E0 + vaccinated exceeds the population in some groups")
  }
  zero <- matrix(0, n_group, 2)

  sched <- build_beta_schedule(beta_schedule, dims)

  pars <- list(
    n_group = n_group,
    n_beta_time = length(sched$time),
    beta = NA_real_,
    sigma = 1 / latent_period,
    gamma = 1 / infectious_period,
    contact = contact,
    sus = vs(susceptibility, susceptibility * (1 - ve_infection)),
    trans = vs(transmissibility, transmissibility * (1 - ve_transmission)),
    p_hosp = vs(p_hosp, p_hosp * (1 - ve_hosp)),
    p_death_comm = vs(p_death_comm, p_death_comm * (1 - ve_death)),
    p_death_hosp = vs(p_death_hosp, p_death_hosp * (1 - ve_death)),
    dur_hosp = vs(dur_hosp, dur_hosp_vacc),
    S_ini = vs(S_unvacc, S_vacc),
    E_ini = vs(E0, 0),
    I_ini = vs(I0, 0),
    R_ini = zero,
    beta_time = sched$time,
    beta_mult_value = sched$value,
    # Observation model, only used when fitting (fit_model() replaces these)
    n_obs = 1L,
    obs_map = matrix(1, 1, n_group),
    nb_size = 10)

  N_group <- as.vector(N)
  if (is.null(beta)) {
    beta <- R0 / r0_per_beta(contact, pars$sus, pars$trans, N_group,
                             infectious_period)
  }
  check_range(beta, "beta")
  pars$beta <- beta

  if (sum(I0 + E0) == 0) {
    warning("No initial infections (I0 and E0 are all zero)")
  }

  ret <- structure(
    list(pars = pars,
         beta = beta,
         dims = dims,
         N = N,
         N_group = N_group,
         C_age = C_age,
         M_region = M_region,
         epi = list(susceptibility = susceptibility,
                    transmissibility = transmissibility,
                    p_hosp = p_hosp, p_death_hosp = p_death_hosp,
                    p_death_comm = p_death_comm, vacc_coverage = vacc_coverage,
                    ve_infection = ve_infection,
                    ve_transmission = ve_transmission,
                    ve_hosp = ve_hosp, ve_death = ve_death),
         args = args),
    class = "seirvax_model")
  ret
}

#' Change parameters of an existing model
#'
#' Rebuilds the model with some arguments of [build_model()] changed. Giving
#' `beta` replaces a previously given `R0` and vice versa.
#'
#' @param model A model from [build_model()].
#' @param ... Arguments of [build_model()] to change.
#' @return A new `seirvax_model`.
#' @export
#' @examples
#' m <- build_model(N = 1e5, R0 = 2, I0 = 10)
#' m2 <- update_model(m, R0 = 1.5)
update_model <- function(model, ...) {
  changes <- list(...)
  if (length(changes) > 0 && (is.null(names(changes)) || any(names(changes) == ""))) {
    stop("All arguments to update_model() must be named")
  }
  unknown <- setdiff(names(changes), names(formals(build_model)))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown argument(s): %s", paste(unknown, collapse = ", ")))
  }
  args <- model$args
  if ("beta" %in% names(changes)) args$R0 <- NULL
  if ("R0" %in% names(changes)) args$beta <- NULL
  args[names(changes)] <- changes
  do.call(build_model, args)
}

#' @export
print.seirvax_model <- function(x, ...) {
  cat("<seirvax_model>\n")
  cat(sprintf("  %d age group(s) x %d region(s); population %s\n",
              length(x$dims$age), length(x$dims$region),
              format(sum(x$N), big.mark = ",", scientific = FALSE)))
  cat(sprintf("  Age groups: %s\n", paste(x$dims$age, collapse = ", ")))
  cat(sprintf("  Regions:    %s\n", paste(x$dims$region, collapse = ", ")))
  cat(sprintf("  R0 = %.3g, beta = %.4g, R_eff(0) = %.3g\n",
              R0(x), x$beta, R_eff(x)))
  cat(sprintf("  Latent period %.3g days, infectious period %.3g days\n",
              1 / x$pars$sigma, 1 / x$pars$gamma))
  cat(sprintf("  Initial infected: %s (I0 + E0); vaccinated: %s\n",
              format(sum(x$pars$I_ini + x$pars$E_ini), big.mark = ","),
              format(sum(x$pars$S_ini[, 2]), big.mark = ",")))
  if (x$pars$n_beta_time > 1) {
    cat(sprintf("  Contact multipliers change at time(s): %s\n",
                paste(x$pars$beta_time[-1], collapse = ", ")))
  }
  cat("  Parameters by age (averaged over regions, population-weighted):\n")
  w <- x$N / rowSums(x$N)
  tab <- vapply(x$epi, function(p) rowSums(p * w), numeric(length(x$dims$age)))
  tab <- as.data.frame(matrix(tab, nrow = length(x$dims$age),
                              dimnames = list(x$dims$age, names(x$epi))))
  tab <- t(tab)
  print(noquote(apply(tab, c(1, 2), function(v) format(signif(v, 3)))))
  invisible(x)
}
