test_that("read_demography returns tidy data with ordered age groups", {
  d <- read_demography()
  expect_setequal(names(d), c("country", "region", "age_group", "population",
                              "share_of_region", "area_km2", "density_per_km2",
                              "pop_total", "profile"))
  expect_setequal(unique(d$country), c("Sweway", "Finmark", "Norland", "Denden"))
  expect_equal(uniqueN(d$region), 20)
  expect_true(is.ordered(d$age_group))
  expect_equal(levels(d$age_group)[1:3], c("0-4", "5-9", "10-14"))
  expect_equal(tail(levels(d$age_group), 1), "95+")
  # population by region matches the region total in the data
  tot <- d[, .(p = sum(population), n = pop_total[1]), by = region]
  expect_equal(tot$p, tot$n)
})

test_that("region_info has one row per region", {
  info <- region_info()
  expect_equal(nrow(info), 20)
  expect_equal(nrow(region_info(country = "Sweway")), 5)
  expect_equal(region_info(regions = c("S2", "D1"))$region, c("S2", "D1"))
  expect_error(region_info(country = "Nowhere"), "Unknown country")
  expect_error(region_info(regions = "X9"), "Unknown region")
})

test_that("population_matrix selects and aggregates", {
  d <- read_demography()
  full <- population_matrix(demography = d)
  expect_equal(dim(full), c(20, 20))
  expect_equal(sum(full), sum(d$population))

  one <- population_matrix("Sweway", age_breaks = c(0, 20, 45, 65))
  expect_equal(rownames(one), c("0-19", "20-44", "45-64", "65+"))
  expect_equal(colnames(one), paste0("S", 1:5))
  expect_equal(sum(one), sum(d[country == "Sweway"]$population))
  # 0-19 is the first four five-year bands
  expect_equal(unname(one["0-19", "S1"]),
               sum(d[region == "S1" & age_group %in% c("0-4", "5-9", "10-14", "15-19")]$population))

  # region order follows the 'regions' argument
  two <- population_matrix(regions = c("D3", "S1"))
  expect_equal(colnames(two), c("D3", "S1"))

  # one age group, countries instead of regions
  by_country <- population_matrix(by_region = FALSE, age_breaks = 0)
  expect_equal(dim(by_country), c(1, 4))
  expect_equal(rownames(by_country), "0+")
  expect_equal(sum(by_country), sum(d$population))

  expect_error(population_matrix(age_breaks = c(0, 22)), "start of an age band")
  expect_error(population_matrix(age_breaks = c(20, 65)), "must be 0")
})

test_that("reading data does not modify the caller's table", {
  d <- read_demography()
  before <- copy(d)
  population_matrix("Sweway", age_breaks = c(0, 65), demography = d)
  region_info(regions = "S3", demography = d)
  expect_equal(d, before)
})

test_that("the population matrix works with build_model", {
  N <- population_matrix("Sweway", age_breaks = c(0, 20, 65))
  m <- build_model(N = N, C_age = matrix(1, 3, 3), M_region = diag(5),
                   R0 = 2, I0 = 1, p_hosp = 0.01)
  expect_equal(m$dims$age, c("0-19", "20-64", "65+"))
  expect_equal(m$dims$region, paste0("S", 1:5))
  expect_equal(sum(m$N), sum(N))
})

test_that("geography helpers read and filter the map layers", {
  skip_if_not_installed("sf")
  regions <- read_geography("regions", country = "Sweway")
  expect_equal(nrow(regions), 5)
  expect_true(all(c("region", "country", "POP_TOTAL") %in% names(regions)))
  expect_equal(nrow(read_geography("settlements", regions = "S1")),
               sum(read_geography("settlements")$region == "S1"))
  expect_gt(nrow(read_geography("cells", country = "Denden")), 100)
  expect_error(read_geography("regions", country = "Nowhere"), "No features")
})

test_that("region_distances is a symmetric km matrix in the given order", {
  skip_if_not_installed("sf")
  d <- region_distances(regions = c("S2", "S1", "D1"))
  expect_equal(dimnames(d), list(c("S2", "S1", "D1"), c("S2", "S1", "D1")))
  expect_equal(unname(diag(d)), rep(0, 3))
  expect_equal(d, t(d))
  expect_true(all(d[upper.tri(d)] > 0) && max(d) < 5000)   # kilometres
})

test_that("plot_region_map accepts vectors and data.frames", {
  skip_if_not_installed("sf")
  expect_s3_class(plot_region_map(country = "Sweway"), "ggplot")
  v <- c(S1 = 1, S2 = 2, S3 = 3, S4 = 4, S5 = 5)
  expect_s3_class(plot_region_map(v), "ggplot")
  expect_s3_class(plot_region_map(data.frame(region = names(v), attack_rate = v)),
                  "ggplot")
  expect_error(plot_region_map(c(1, 2)), "named with regions")
})
