# A small two-age, two-region model used across tests.
small_model <- function(...) {
  N <- matrix(c(20000, 30000, 10000, 15000), 2,
              dimnames = list(c("young", "old"), c("A", "B")))
  C_age <- matrix(c(3, 1, 1, 2), 2)
  M_region <- matrix(c(0.9, 0.1, 0.1, 0.9), 2)
  args <- list(N = N, C_age = C_age, M_region = M_region, R0 = 2,
               I0 = matrix(c(10, 0, 0, 0), 2), p_hosp = c(0.01, 0.1),
               p_death_hosp = c(0.02, 0.2), p_death_comm = c(0, 0.01))
  args[names(list(...))] <- list(...)
  do.call(build_model, args)
}

total_alive_and_dead <- function(res) {
  d <- get_output(res, c("S", "E", "I", "H", "R", "D"))
  d[, .(total = sum(value)), by = .(sim, time)]$total
}
