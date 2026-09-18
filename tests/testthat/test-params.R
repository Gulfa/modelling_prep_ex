dims <- list(age = c("a1", "a2", "a3"), region = c("r1", "r2"))

test_that("expand_param recycles scalars, age vectors, region vectors and matrices", {
  expect_equal(unname(expand_param(2, "x", dims)), matrix(2, 3, 2))
  expect_equal(unname(expand_param(1:3, "x", dims)), matrix(1:3, 3, 2))
  expect_equal(unname(expand_param(c(5, 6), "x", dims)),
               matrix(c(5, 6), 3, 2, byrow = TRUE))
  expect_equal(unname(expand_param(c(r2 = 6, r1 = 5), "x", dims)),
               matrix(c(5, 6), 3, 2, byrow = TRUE))
  m <- matrix(1:6, 3, 2)
  expect_equal(unname(expand_param(m, "x", dims)), m)
  expect_error(expand_param(1:4, "x", dims), "vector by age")
  expect_error(expand_param(matrix(1, 2, 3), "x", dims), "3 x 2")
})

test_that("as.vector of an age x region matrix is in backend group order", {
  m <- expand_param(matrix(1:6, 3, 2), "x", dims)
  # group (r - 1) * n_age + a
  expect_equal(as.vector(m)[(2 - 1) * 3 + 3], m["a3", "r2"])
})

test_that("combine_mixing is the Kronecker product with age fastest", {
  C_age <- matrix(c(1, 2, 3, 4), 2)
  M <- matrix(c(0.7, 0.2, 0.3, 0.8), 2)
  C <- combine_mixing(C_age, M)
  # age b = 2 in region s = 1 seen from age a = 1 in region r = 2
  expect_equal(C[(2 - 1) * 2 + 1, (1 - 1) * 2 + 2], C_age[1, 2] * M[2, 1])
})

test_that("build_model validates its inputs", {
  expect_error(build_model(N = 100, I0 = 1), "exactly one")
  expect_error(build_model(N = 100, R0 = 2, beta = 1, I0 = 1), "exactly one")
  expect_error(build_model(N = c(a = 10, b = 10), R0 = 2, I0 = 1), "C_age")
  expect_error(build_model(N = 100, R0 = 2, I0 = 1, p_hosp = 1.5), "p_hosp")
  expect_error(build_model(N = 100, R0 = 2, I0 = 200), "exceeds")
  expect_error(small_model(M_region = matrix(0.4, 2, 2)), "sum to 1")
  expect_warning(build_model(N = 100, R0 = 2), "No initial infections")
})

test_that("vaccinated parameters are derived from VE", {
  m <- small_model(vacc_coverage = 0.5, ve_infection = 0.6, ve_hosp = 0.5,
                   ve_transmission = 0.2, ve_death = 0.25, dur_hosp_vacc = 4)
  p <- m$pars
  expect_equal(p$sus[, 2], p$sus[, 1] * 0.4)
  expect_equal(p$trans[, 2], p$trans[, 1] * 0.8)
  expect_equal(p$p_hosp[, 2], p$p_hosp[, 1] * 0.5)
  expect_equal(p$p_death_hosp[, 2], p$p_death_hosp[, 1] * 0.75)
  expect_equal(p$dur_hosp[, 2], rep(4, 4))
  expect_equal(rowSums(p$S_ini + p$I_ini), as.vector(m$N))
})

test_that("update_model changes arguments and swaps R0/beta", {
  m <- small_model()
  m2 <- update_model(m, ve_infection = 0.3, vacc_coverage = 0.5)
  expect_equal(m2$pars$sus[, 2], m2$pars$sus[, 1] * 0.7)
  m3 <- update_model(m, beta = 0.1)
  expect_equal(m3$beta, 0.1)
  m4 <- update_model(m3, R0 = 3)
  expect_equal(R0(m4), 3)
  expect_error(update_model(m, not_an_arg = 1), "Unknown")
})
