test_that("get_output aggregates consistently", {
  m <- small_model(vacc_coverage = 0.3)
  res <- run_model(m, t_max = 50)
  full <- get_output(res, "I", by = c("age", "region", "vacc_status"))
  tot <- get_output(res, "I")
  expect_equal(full[, .(v = sum(value)), by = time]$v, tot$value)
  expect_equal(levels(get_output(res, "I", by = "region")$region), c("A", "B"))
  expect_error(get_output(res, "nope"), "Unknown variable")
  expect_error(get_output(res, "I", by = "county"), "subset")
})

test_that("summary and plot work for both model types", {
  m <- small_model()
  det <- run_model(m, t_max = 100)
  sto <- run_model(m, "stochastic", t_max = 100, n_sim = 4, seed = 1)
  s1 <- summary(det, by = "age")
  expect_setequal(names(s1), c("age", "measure", "value"))
  s2 <- summary(sto)
  expect_setequal(names(s2), c("measure", "mean", "median", "lower", "upper"))
  expect_s3_class(plot(det, "H", by = c("age", "region")), "ggplot")
  expect_s3_class(plot(sto, "inc_hosp", by = "age", n_traces = 2), "ggplot")
})

test_that("run_model can keep a subset of variables", {
  res <- run_model(small_model(), t_max = 20, vars = c("H", "inc_deaths"))
  expect_equal(names(res$states), c("H", "inc_deaths"))
  expect_error(run_model(small_model(), t_max = 20, vars = "nope"), "Unknown")
})
