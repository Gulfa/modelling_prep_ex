small_line_list <- function(...) {
  args <- list(country = "Sweway", R0 = 2, t_max = 50, seed = 1, min_cases = 1,
               max_tries = 1)
  args[names(list(...))] <- list(...)
  do.call(simulate_line_list, args)
}

test_that("the line list has consistent records and dates", {
  ll <- small_line_list()
  expect_named(ll, c("id", "date_onset", "age_group", "region", "vaccinated",
                     "hospitalised", "date_hospital", "died", "date_death"))
  expect_gt(nrow(ll), 0)
  expect_equal(ll$id, seq_len(nrow(ll)))
  expect_false(is.unsorted(ll$date_onset))

  # dates line up with the flags
  expect_true(all(!is.na(ll$date_hospital) == ll$hospitalised))
  expect_true(all(!is.na(ll$date_death) == ll$died))
  expect_true(all(ll$date_onset <= ll$date_hospital, na.rm = TRUE))
  expect_true(all(ll$date_hospital <= ll$date_death, na.rm = TRUE))
  expect_true(all(ll$date_onset <= ll$date_death, na.rm = TRUE))
  expect_s3_class(ll$date_onset, "Date")

  # only regions of the country, and the model's age bands
  expect_true(all(ll$region %in% paste0("S", 1:5)))
  expect_setequal(levels(ll$age_group), rownames(population_matrix("Sweway")))
})

test_that("severe cases are always detected and p_detect controls the rest", {
  severe_only <- small_line_list(p_detect = 0)
  expect_true(all(severe_only$hospitalised | severe_only$died))

  few <- small_line_list(p_detect = 0.05)
  many <- small_line_list(p_detect = 0.5)
  expect_gt(nrow(many), nrow(few))
  # the severe cases are the same in both (same seed, same epidemic)
  expect_equal(sum(few$hospitalised), sum(many$hospitalised))
})

test_that("age_breaks only changes the reporting", {
  ll <- small_line_list(age_breaks = c(0, 20, 65))
  expect_equal(levels(ll$age_group), c("0-19", "20-64", "65+"))
  expect_equal(nrow(ll), nrow(small_line_list()))
})

test_that("R0, seeding and the truth attribute", {
  ll <- small_line_list(R0 = NULL, R0_range = c(1.5, 1.6))
  truth <- attr(ll, "truth")
  expect_true(truth$R0 >= 1.5 && truth$R0 <= 1.6)
  expect_equal(R0(truth$model), truth$R0)
  expect_equal(truth$country, "Sweway")

  seeded <- small_line_list(seed_region = "S3")
  expect_equal(attr(seeded, "truth")$model$pars$I_ini[, 1][
    rep(colnames(population_matrix("Sweway")), each = 20) != "S3"], rep(0, 80))
  expect_error(small_line_list(seed_region = "D1"), "must be a region")
})

test_that("results are reproducible and start_date can be dropped", {
  expect_equal(small_line_list(), small_line_list())
  days <- small_line_list(start_date = NULL)
  expect_type(days$date_onset, "double")
  expect_gte(min(days$date_onset), 0)
})

test_that("a faded-out epidemic warns rather than failing", {
  expect_warning(simulate_line_list("Sweway", R0 = 0.5, t_max = 20, n_seed = 1,
                                    seed = 2, min_cases = 1000, max_tries = 2),
                 "faded out")
})

test_that("admissions_from_line_list counts every admission", {
  ll <- small_line_list()
  d <- admissions_from_line_list(ll)
  expect_named(d, c("time", "admissions"))
  expect_equal(sum(d$admissions), sum(ll$hospitalised))
  expect_equal(d$time, seq_len(max(d$time)))     # includes days with no cases
  expect_equal(min(d$time), 1)

  by_region <- admissions_from_line_list(ll, by = "region")
  expect_equal(sum(by_region$admissions), sum(ll$hospitalised))
  expect_setequal(names(by_region), c("time", "region", "admissions"))

  # an earlier start date shifts the series
  start <- min(ll$date_hospital, na.rm = TRUE) - 10
  shifted <- admissions_from_line_list(ll, start_date = start)
  expect_equal(max(shifted$time), max(d$time) + 10)
  expect_equal(sum(shifted$admissions), sum(d$admissions))

  expect_error(admissions_from_line_list(ll, by = "county"), "subset")
  expect_error(admissions_from_line_list(ll, start_date = max(ll$date_hospital,
                                                              na.rm = TRUE)),
               "after the first admission")
  expect_error(admissions_from_line_list(ll, start_date = NA), "non-missing")
  expect_error(admissions_from_line_list(ll[0]), "no hospital admissions")
})

test_that("the generated admissions can be fitted", {
  truth_R0 <- 2
  ll <- small_line_list(R0 = truth_R0, t_max = 60, p_detect = 0)
  d <- admissions_from_line_list(ll, start_date = min(ll$date_onset))
  truth <- attr(ll, "truth")$model
  # a coarser model of the same epidemic: same natural history and seeding as
  # the generator, so that R0 is the main thing left to estimate
  N <- population_matrix("Sweway", age_breaks = c(0, 20, 65))
  seeds <- matrix(rowSums(truth$pars$I_ini), length(truth$dims$age),
                  dimnames = list(truth$dims$age, truth$dims$region))
  seeded <- colnames(seeds)[colSums(seeds) > 0]
  I0 <- matrix(0, 3, 5, dimnames = dimnames(N))
  I0["20-64", seeded] <- sum(seeds)
  expect_equal(sum(I0), sum(truth$pars$I_ini))
  m <- build_model(N = N, C_age = age_contact_matrix(rownames(N)),
                   M_region = truth$M_region, R0 = 3,
                   latent_period = 1 / truth$pars$sigma,
                   infectious_period = 1 / truth$pars$gamma,
                   p_hosp = c(0.005, 0.02, 0.1), I0 = I0)
  fit <- fit_model(m, d, pars = "R0")
  # structure is coarser than the generator's, so allow a wide band
  expect_gt(fit$estimates$estimate, 0.6 * truth_R0)
  expect_lt(fit$estimates$estimate, 1.6 * truth_R0)
})

test_that("age_contact_matrix is symmetric with the right labels", {
  ages <- c("0-19", "20-44", "45-64", "65+")
  m <- age_contact_matrix(ages)
  expect_equal(dimnames(m), list(ages, ages))
  expect_equal(m, t(m))
  expect_true(all(diag(m) == apply(m, 1, max)))   # most contacts within age
})

test_that("the epidemiological parameters are arguments", {
  ll <- small_line_list(latent_period = 3, infectious_period = 5,
                        dur_hosp = 12, ve_infection = 0.8)
  p <- attr(ll, "truth")$model$pars
  expect_equal(1 / p$sigma, 3)
  expect_equal(1 / p$gamma, 5)
  expect_equal(unique(p$dur_hosp[, 1]), 12)
  expect_equal(p$sus[, 2], p$sus[, 1] * 0.2)
})

test_that("age-varying parameters take numbers, vectors or functions", {
  n_bands <- nrow(population_matrix("Sweway"))

  flat <- small_line_list(p_hosp = 0.02)
  expect_equal(unique(attr(flat, "truth")$model$pars$p_hosp[, 1]), 0.02)

  by_band <- small_line_list(susceptibility = seq(0.5, 1, length.out = n_bands))
  sus <- attr(by_band, "truth")$model$pars$sus[, 1]
  expect_equal(length(unique(sus)), n_bands)

  # a function is evaluated at the midpoint of each band
  f <- small_line_list(vacc_coverage = function(age) ifelse(age >= 65, 0.8, 0))
  m <- attr(f, "truth")$model
  cov <- m$epi$vacc_coverage[, 1]
  expect_equal(unname(cov[age_lower(names(cov)) >= 65]), rep(0.8, 7))
  expect_equal(unname(cov[age_lower(names(cov)) < 60]), rep(0, 12))
  expect_gt(sum(m$pars$S_ini[, 2]), 0)
  expect_true(any(f$vaccinated))

  expect_error(small_line_list(p_hosp = c(0.1, 0.2)), "one value per")
  expect_error(small_line_list(dur_hosp = "long"), "must be a number")
})

test_that("no vaccination by default", {
  m <- attr(small_line_list(), "truth")$model
  expect_equal(sum(m$pars$S_ini[, 2]), 0)
  expect_equal(R_eff(m), R0(m), tolerance = 1e-3)
})

test_that("seeding is controlled by n_seed and seed_ages", {
  ll <- small_line_list(n_seed = 50, seed_ages = c(60, 100))
  m <- attr(ll, "truth")$model
  seeds <- matrix(rowSums(m$pars$I_ini), length(m$dims$age),
                  dimnames = list(m$dims$age, m$dims$region))
  expect_equal(sum(seeds), 50)
  expect_equal(sum(seeds[age_lower(rownames(seeds)) < 60, ]), 0)
  expect_error(small_line_list(seed_ages = c(60, 20)), "youngest first")
  expect_error(small_line_list(seed_ages = c(100, 120)), "No age band")
})
