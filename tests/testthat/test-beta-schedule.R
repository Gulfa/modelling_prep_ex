dims <- list(age = c("a1", "a2"), region = c("r1", "r2"))

test_that("no schedule gives a constant multiplier of 1", {
  s <- build_beta_schedule(NULL, dims)
  expect_equal(s$time, 0)
  expect_equal(s$value, matrix(1, 4, 1))
})

test_that("a global schedule applies to every group", {
  s <- build_beta_schedule(data.frame(time = c(10, 20), mult = c(0.5, 1)), dims)
  expect_equal(s$time, c(0, 10, 20))
  expect_equal(s$value, matrix(rep(c(1, 0.5, 1), each = 4), 4))
})

test_that("region and age specific rows override only matching groups", {
  sched <- data.frame(time = c(10, 20, 30), mult = c(0.5, 0.2, 1),
                      region = c(NA, "r2", NA), age = c(NA, "a1", NA))
  s <- build_beta_schedule(sched, dims)
  # groups: a1r1, a2r1, a1r2, a2r2
  expect_equal(s$value[, 3], c(0.5, 0.5, 0.2, 0.5))
  expect_equal(s$value[, 4], rep(1, 4))
  expect_error(build_beta_schedule(data.frame(time = 1, mult = 1, region = "x"), dims),
               "Unknown region")
})

test_that("the multiplier changes transmission in the model", {
  m <- small_model(beta_schedule = data.frame(time = 20, mult = 0))
  res <- run_model(m, t_max = 60)
  inc <- get_output(res, "inc_infections")
  # no new infections once E has emptied is too strict; infections stop at t > 20
  expect_true(all(inc[time > 20]$value < 1e-8))
  expect_true(sum(inc[time <= 20]$value) > 0)
})

test_that("regional lockdown with isolated regions only affects that region", {
  m <- small_model(M_region = diag(2), I0 = matrix(10, 2, 2),
                   beta_schedule = data.frame(time = 0, mult = 0, region = "B"))
  res <- run_model(m, t_max = 50)
  cum <- get_output(res, "cum_infections", by = "region")[time == 50]
  expect_equal(cum[region == "B"]$value, 0)
  expect_gt(cum[region == "A"]$value, 100)
})

test_that("a targeted multiplier scales that group's contacts with everyone", {
  # "old" at 0.5 should equal scaling C_age by min(m_a, m_b) with the same beta
  m <- small_model(beta_schedule = data.frame(time = 0, mult = 0.5, age = "old"))
  C_age <- matrix(c(3, 1, 1, 2), 2)
  C_scaled <- C_age * outer(c(1, 0.5), c(1, 0.5), pmin)   # old-young halved too
  ref <- small_model(C_age = C_scaled, R0 = NULL, beta = m$beta)

  expect_equal(R_eff(m), R_eff(ref), tolerance = 1e-10)
  a <- get_output(run_model(m, t_max = 150), c("H", "cum_infections"),
                  by = c("age", "region"))
  b <- get_output(run_model(ref, t_max = 150), c("H", "cum_infections"),
                  by = c("age", "region"))
  expect_equal(a$value, b$value, tolerance = 1e-6)
})

test_that("a global multiplier scales R_eff by exactly that factor", {
  m <- small_model(beta_schedule = data.frame(time = 10, mult = 0.6))
  expect_equal(R_eff(m, time = 10), 0.6 * R_eff(m, time = 0))
})
