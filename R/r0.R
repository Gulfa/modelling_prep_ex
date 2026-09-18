# Next-generation matrix (NGM) calculations.
#
# With the force of infection used in the models, one infected person in
# group j, stratum w, infects on average
#
#   K[(i,v), (j,w)] = beta * min(m_i, m_j) * sus[i,v] * C[i,j] * trans[j,w]
#                     * S[i,v] / N_j * infectious_period
#
# people in group i, stratum v. Strata are stacked unvaccinated first.
ngm <- function(beta, contact, sus, trans, S, N_group, mult, infectious_period) {
  a <- as.vector(sus * S)
  b <- as.vector(trans / N_group)
  C2 <- kronecker(matrix(1, 2, 2), contact * outer(mult, mult, pmin))
  beta * infectious_period * (a * C2) * rep(b, each = length(a))
}

spectral_radius <- function(K) {
  max(Mod(eigen(K, only.values = TRUE)$values))
}

# R0 per unit beta: fully susceptible, unvaccinated population, no contact
# multipliers.
r0_per_beta <- function(contact, sus, trans, N_group, infectious_period) {
  S <- cbind(N_group, 0)
  spectral_radius(ngm(1, contact, sus, trans, S, N_group,
                      rep(1, length(N_group)), infectious_period))
}

#' Reproduction numbers of a model
#'
#' `R0` is the basic reproduction number: the spectral radius of the
#' next-generation matrix for a fully susceptible, unvaccinated population
#' without contact multipliers. `R_eff` uses the susceptible population
#' (including vaccinated) and the contact multipliers at time `time`,
#' assuming nobody has been infected yet.
#'
#' Note that the NGM R0 gives the right *generation-to-generation* growth once
#' the epidemic has settled into the dominant eigenvector's age/region
#' distribution. If the initial infections are distributed differently (for
#' example seeded in one region), early growth can be faster or slower than
#' R0 suggests.
#'
#' @param model A model from [build_model()].
#' @param time Time at which to evaluate the contact multipliers (R_eff only).
#' @return A single number.
#' @export
R0 <- function(model) {
  p <- model$pars
  model$beta * r0_per_beta(p$contact, p$sus, p$trans, model$N_group,
                           1 / p$gamma)
}

#' @rdname R0
#' @export
R_eff <- function(model, time = 0) {
  p <- model$pars
  k <- findInterval(time, p$beta_time)
  S <- p$S_ini + p$E_ini + p$I_ini
  spectral_radius(ngm(p$beta, p$contact, p$sus, p$trans, S, model$N_group,
                      p$beta_mult_value[, k], 1 / p$gamma))
}
