# Noise-free admissions from the deterministic model
det_admissions <- function(model, t_max, by = character(0)) {
  d <- get_output(run_model(model, t_max = t_max), "inc_hosp", by = by)[time > 0]
  d[, admissions := round(value)]
  d[, c("time", by, "admissions"), with = FALSE]
}

test_that("prepare_fit_data builds series, obs_map and dust2 data", {
  m <- small_model()
  d <- data.frame(time = c(1, 1, 2, 3), region = c("A", "B", "A", "B"),
                  admissions = c(0, 1, 2, 3))
  obs <- prepare_fit_data(d, m)
  expect_equal(obs$series$label, c("A", "B"))
  # groups: young A, old A, young B, old B
  expect_equal(obs$obs_map, rbind(c(1, 1, 0, 0), c(0, 0, 1, 1)))
  expect_equal(obs$data_dust$time, 1:3)
  expect_equal(obs$data_dust$admissions[[2]], c(2, NA))

  national <- prepare_fit_data(data.frame(time = 1:3, admissions = 0), m)
  expect_equal(national$obs_map, matrix(1, 1, 4))

  expect_error(prepare_fit_data(data.frame(time = 0, admissions = 1), m), ">= 1")
  expect_error(prepare_fit_data(data.frame(time = 1, admissions = 1.5), m), "whole")
  expect_error(prepare_fit_data(data.frame(time = 1, admissions = 1, region = "X"), m),
               "Unknown region")
  expect_error(prepare_fit_data(data.frame(time = c(1, 1), admissions = 1), m),
               "more than one row")
})

test_that("optim recovers R0 from noise-free data", {
  truth <- small_model(R0 = 1.7)
  d <- det_admissions(truth, 150)
  fit <- fit_model(update_model(truth, R0 = 2.5), d)
  expect_equal(fit$estimates$estimate, 1.7, tolerance = 0.01)
  expect_equal(R0(fit$model), fit$estimates$estimate)
  expect_true(fit$estimates$lower < 1.7 && fit$estimates$upper > 1.7)
})

test_that("optim recovers R0 and hosp_scale from regional data", {
  truth <- small_model(R0 = 1.7)
  d <- det_admissions(update_model(truth, p_hosp = c(0.02, 0.2)), 150,
                      by = "region")
  fit <- fit_model(truth, d, pars = c("R0", "hosp_scale"))
  expect_equal(fit$estimates$estimate, c(1.7, 2), tolerance = 0.02)
  expect_equal(fit$model$pars$p_hosp[, 1], rep(c(0.01, 0.1), 2) * fit$estimates$estimate[2])
})

test_that("custom priors and starting values are used", {
  m <- small_model()
  d <- det_admissions(m, 60)
  fit <- fit_model(m, d, prior = c(R0 = "Uniform(1.9, 3)"), initial = c(R0 = 2.5))
  expect_gte(fit$estimates$estimate, 1.9)   # truth (2) is near the prior edge
  expect_error(fit_model(m, d, prior = c(beta = "Uniform(0, 1)")), "named")
  expect_error(fit_model(m, d, pars = "gamma"), "one or more of")
  expect_error(fit_model(m, d, initial = c(R0 = 10)), "zero prior")
})

test_that("mcmc returns samples around the truth", {
  truth <- small_model(R0 = 1.7)
  set.seed(1)
  d <- det_admissions(truth, 120)
  d[, admissions := rnbinom(.N, size = 10, mu = admissions + 0.1)]
  fit <- suppressMessages(fit_model(truth, d, pars = c("R0", "nb_size"),
                                    method = "mcmc", n_steps = 400,
                                    n_chains = 2, seed = 1))
  expect_equal(dim(fit$samples$pars), c(2, 300, 2))
  est <- fit$estimates
  expect_true(est$lower[1] < 1.7 && est$upper[1] > 1.7)
  expect_gt(fit$acceptance, 0.05)
  expect_equal(fit$nb_size, est$estimate[2])
})

test_that("pmcmc runs on the stochastic model", {
  truth <- build_model(N = 20000, R0 = 2, I0 = 20, p_hosp = 0.05)
  d <- det_admissions(truth, 40)
  fit <- suppressMessages(fit_model(truth, d, method = "pmcmc", n_steps = 20,
                                    n_chains = 1, n_particles = 20, dt = 0.25,
                                    seed = 1))
  expect_equal(dim(fit$samples$pars), c(1, 15, 1))
  expect_true(all(is.finite(fit$samples$density)))
})

test_that("project_fit attaches the draws it used and can be seeded", {
  truth <- small_model(R0 = 1.7)
  set.seed(1)
  d <- det_admissions(truth, 120)
  d[, admissions := rnbinom(.N, size = 10, mu = admissions + 0.1)]
  fit <- suppressMessages(fit_model(truth, d, pars = c("R0", "nb_size"),
                                    method = "mcmc", n_steps = 400,
                                    n_chains = 2, seed = 1))
  p1 <- project_fit(fit, from = 0, t_max = 60, n_draws = 10, seed = 42)
  draws <- p1$draws
  expect_equal(dim(draws), c(10, 2))
  expect_equal(colnames(draws), c("R0", "nb_size"))
  expect_equal(p1$n_sim, 10)              # one simulation per draw
  expect_equal(p1$draw, 1:10)

  # the same seed gives the same draws and the same trajectories
  p2 <- project_fit(fit, from = 0, t_max = 60, n_draws = 10, seed = 42)
  expect_equal(p2$draws, draws)
  expect_equal(get_output(p2, "inc_hosp"), get_output(p1, "inc_hosp"))
  p3 <- project_fit(fit, from = 0, t_max = 60, n_draws = 10, seed = 7)
  expect_false(isTRUE(all.equal(p3$draws, draws)))

  # each trajectory really comes from its own draw
  one <- project_fit(fit, from = 0, t_max = 60, n_draws = 1, seed = 3)
  m <- fit_rebuild(fit$base_model, one$draws[1, ])
  direct <- get_output(run_model(m, t_max = 60), "inc_hosp")
  expect_equal(get_output(one, "inc_hosp")$value, direct$value)

  # plot() uses those same draws rather than drawing a second set
  expect_s3_class(plot(fit, n_draws = 10, seed = 42), "ggplot")
  a <- plot(fit, n_draws = 10, seed = 1)$layers[[1]]$data
  b <- plot(fit, n_draws = 10, seed = 1)$layers[[1]]$data
  expect_equal(a, b)
})

test_that("projections can use a scenario model", {
  truth <- small_model(R0 = 1.7)
  d <- det_admissions(truth, 100)
  fit <- fit_model(truth, d, pars = "R0")

  base <- project_fit(fit, t_max = 200, n_draws = 5, vars = "cum_hosp", seed = 1)
  stopped <- project_fit(fit, model = update_model(fit$base_model,
                           beta_schedule = data.frame(time = 50, mult = 0)),
                         t_max = 200, n_draws = 5, vars = "cum_hosp", seed = 1)
  # same draws, different model
  expect_equal(stopped$draws, base$draws)
  expect_lt(max(get_output(stopped, "cum_hosp")$value),
            max(get_output(base, "cum_hosp")$value))

  # the drawn parameters are still applied to the scenario model
  one <- project_fit(fit, model = update_model(fit$base_model, p_hosp = 0.5),
                     from = 0, t_max = 50, n_draws = 1, vars = "cum_hosp",
                     seed = 2)
  direct <- get_output(
    run_model(update_model(fit$base_model, p_hosp = 0.5,
                           R0 = one$draws[1, "R0"]), t_max = 50),
    "cum_hosp")
  expect_equal(get_output(one, "cum_hosp")$value, direct$value)

  # project_summary takes it too
  s <- project_summary(fit, model = update_model(fit$base_model,
                         beta_schedule = data.frame(time = 50, mult = 0)),
                       t_max = 200, n_draws = 5, measures = "cum_hosp", seed = 1)
  expect_lt(s$median, project_summary(fit, t_max = 200, n_draws = 5,
                                      measures = "cum_hosp", seed = 1)$median)

  expect_error(project_fit(fit, model = "nope"), "build_model")
  expect_error(project_fit(fit, model = build_model(N = 1e4, R0 = 2, I0 = 5)),
               "same age groups and regions")
})

test_that("projections, plots and printing work", {
  truth <- small_model(R0 = 1.7)
  d <- det_admissions(truth, 60, by = "region")
  fit <- fit_model(truth, d)
  p <- project_fit(fit, t_max = 100, n_draws = 5)
  expect_equal(p$draw, 1:5)
  expect_equal(max(p$time), 100)
  expect_setequal(names(get_output(p, "H", by = "region")),
                  c("time", "sim", "region", "variable", "value"))
  expect_s3_class(plot(fit, n_draws = 5), "ggplot")
  expect_error(plot(fit, "trace"), "No MCMC samples")
  expect_output(print(fit), "maximum a posteriori")
  expect_equal(summary(fit), fit$estimates)
  expect_setequal(fit_parameters()$name, c("R0", "hosp_scale", "nb_size"))
})

test_that("a projection behaves like a model run", {
  truth <- small_model(R0 = 1.7)
  fit <- fit_model(truth, det_admissions(truth, 80), pars = "R0")
  p <- project_fit(fit, t_max = 150, n_draws = 6, seed = 1)
  expect_s3_class(p, "seirvax_projection")
  expect_s3_class(p, "seirvax_run")

  expect_s3_class(plot(p, var = "H"), "ggplot")
  expect_s3_class(plot(p, var = "inc_hosp", by = "region"), "ggplot")
  expect_output(print(p), "seirvax_projection")

  s <- summary(p)
  expect_setequal(names(s), c("measure", "mean", "median", "lower", "upper"))
  expect_setequal(s$measure, c("cum_infections", "cum_hosp", "D",
                               "peak_hosp_occupancy", "peak_hosp_day"))
  # summarised over draws, so there is spread
  expect_gt(s[measure == "cum_hosp"]$upper, s[measure == "cum_hosp"]$lower)

  # several simulations per draw stack up as simulations
  pm <- project_fit(fit, t_max = 120, n_draws = 4, type = "stochastic",
                    n_sim = 3, seed = 1)
  expect_equal(pm$n_sim, 12)
  expect_equal(pm$draw, rep(1:4, each = 3))
})

test_that("project_summary reports outcomes by age and region", {
  truth <- small_model(R0 = 1.7)
  d <- det_admissions(truth, 100)
  fit <- fit_model(truth, d, pars = "R0")

  nat <- project_summary(fit, t_max = 200, n_draws = 10, seed = 1)
  expect_setequal(names(nat), c("measure", "mean", "median", "lower", "upper"))
  expect_setequal(nat$measure, c("cum_infections", "cum_hosp", "D",
                                 "peak_hosp_occupancy", "peak_hosp_day"))
  expect_true(all(nat$lower <= nat$median & nat$median <= nat$upper))

  by_age <- project_summary(fit, by = "age", t_max = 200, n_draws = 10,
                            measures = c("cum_infections", "cum_hosp"),
                            seed = 1)
  expect_setequal(by_age$age, c("young", "old"))
  expect_setequal(by_age$measure, c("cum_infections", "cum_hosp"))
  expect_true(all(c("share_median", "share_lower", "share_upper") %in% names(by_age)))

  # strata add up to the national total (means are additive)
  totals <- by_age[, .(total = sum(mean)), by = measure]
  expect_equal(totals$total,
               nat[match(totals$measure, measure)]$mean, tolerance = 1e-6)
  # shares add to 1 within each measure
  expect_equal(by_age[, sum(share_median), by = measure]$V1, c(1, 1),
               tolerance = 1e-3)

  two <- project_summary(fit, by = c("age", "region"), t_max = 120,
                         n_draws = 5, measures = "D", seed = 1)
  expect_equal(nrow(two), 4)
  expect_equal(sum(two$share_mean), 1, tolerance = 1e-6)
})

test_that("project_summary validates and is reproducible", {
  fit <- fit_model(small_model(R0 = 1.7), det_admissions(small_model(R0 = 1.7), 80),
                   pars = "R0")
  expect_error(project_summary(fit, measures = "cases"), "subset of")
  expect_error(project_summary(fit, by = "county"), "subset")
  a <- project_summary(fit, by = "age", t_max = 100, n_draws = 8, seed = 5)
  b <- project_summary(fit, by = "age", t_max = 100, n_draws = 8, seed = 5)
  expect_equal(a, b)
  # peaks have no meaningful share across strata
  expect_true(all(is.na(a[measure %in% c("peak_hosp_day",
                                         "peak_hosp_occupancy")]$share_mean)))
})

test_that("distribution plots", {
  fit <- fit_model(small_model(R0 = 1.7), det_admissions(small_model(R0 = 1.7), 80),
                   pars = "R0")
  expect_s3_class(plot(fit, "distribution", by = "age", n_draws = 5, seed = 1),
                  "ggplot")
  expect_s3_class(plot(fit, "distribution", by = c("age", "region"),
                       measure = "D", n_draws = 5, seed = 1), "ggplot")
  expect_error(plot(fit, "distribution", by = character(0)), "is needed")

  # observed counts can be overlaid without inherit.aes = FALSE, whatever
  # geom they use
  obs <- data.frame(age = c("young", "old"), observed = c(10, 20))
  base <- plot(fit, "distribution", by = "age", n_draws = 5, seed = 1)
  expect_s3_class(ggplot2::ggplot_build(
    base + ggplot2::geom_point(ggplot2::aes(age, observed), data = obs)),
    "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(
    base + ggplot2::geom_errorbar(ggplot2::aes(x = age, ymin = observed,
                                               ymax = observed), data = obs)),
    "ggplot_built")
  two <- plot(fit, "distribution", by = c("age", "region"), n_draws = 5, seed = 1)
  expect_s3_class(ggplot2::ggplot_build(
    two + ggplot2::geom_point(ggplot2::aes(age, observed), data = obs)),
    "ggplot_built")
  expect_error(plot(fit, "distribution", measure = c("D", "cum_hosp")),
               "single outcome")
})

test_that("projections start after the fit by default", {
  truth <- small_model(R0 = 1.7)
  d <- det_admissions(truth, 80)
  fit <- fit_model(truth, d, pars = "R0")

  p <- project_fit(fit, n_draws = 3, seed = 1)
  expect_equal(min(p$time), max(d$time))          # starts where the data ends
  expect_equal(max(p$time), max(d$time) + 90)     # 90 days by default

  late <- project_fit(fit, t_max = 150, n_draws = 3, seed = 1)
  expect_equal(range(late$time), c(max(d$time), 150))

  # from = 0 gives the whole epidemic, and the overlap is the same trajectory
  full <- project_fit(fit, from = 0, t_max = 150, n_draws = 3, seed = 1)
  expect_equal(min(full$time), 0)
  o <- get_output(full, "inc_hosp")[time >= max(d$time)]
  expect_equal(o$value, get_output(late, "inc_hosp")$value)
  # the epidemic continues rather than restarting at `from`
  expect_gt(min(get_output(late, "cum_hosp")$value), 0)

  expect_error(project_fit(fit, t_max = 10), "must be after")
  expect_error(project_fit(fit, from = -1), "'from' must be")
})

test_that("pmcmc projections continue from the filtered state", {
  truth <- small_model(R0 = 1.8, I0 = matrix(c(20, 0, 0, 0), 2))
  set.seed(2)
  sim <- run_model(truth, "stochastic", t_max = 60, n_sim = 1, seed = 11)
  d <- get_output(sim, "inc_hosp")[time > 0, .(time, admissions = value)]
  fit <- suppressMessages(fit_model(truth, d, pars = "R0", method = "pmcmc",
                                    n_steps = 60, n_chains = 1,
                                    n_particles = 50, dt = 0.25, seed = 1))
  conditioned <- project_fit(fit, t_max = 90, n_draws = 10, type = "stochastic",
                             dt = 0.25, seed = 1)
  expect_equal(min(conditioned$time), max(d$time))

  # starting from the filter's state keeps the handover close to the data,
  # unlike an unconditioned run from time 0
  observed <- d[time == max(time)]$admissions
  handover <- get_output(conditioned, "inc_hosp")[time == max(d$time)]$value
  free <- get_output(project_fit(fit, from = 0, t_max = 90, n_draws = 10,
                                 type = "stochastic", dt = 0.25, seed = 1),
                     "inc_hosp")[time == max(d$time)]$value
  expect_lt(stats::median(abs(handover - observed)),
            stats::median(abs(free - observed)))

  # the filtered state carries the epidemic so far: cumulative counts continue
  expect_gt(min(get_output(conditioned, "cum_hosp")$value), 0)
})
