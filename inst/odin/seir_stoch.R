## Stochastic SEIR-HD model with a parallel vaccinated stratum (tau-leaping).
##
## Discrete-time version of seir_ode.R: the same parameters and states, with
## binomial draws over a step of `dt` days (set when the system is created).
## Outflows from I and H are split with chained binomials, so the competing
## outcomes of a compartment never exceed the number of people leaving it.
##
## The parameter interface and the state declarations must stay identical to
## seir_ode.R (see docs/modifying_the_model.md).

## ---- Parameters -----------------------------------------------------------
n_group     <- parameter(type = "integer", constant = TRUE)
n_beta_time <- parameter(type = "integer", constant = TRUE)

beta  <- parameter()   # transmission rate (calibrated from R0 in R)
sigma <- parameter()   # 1 / latent period
gamma <- parameter()   # 1 / infectious period

# contact[i, j]: contacts per day a person in group i has with group j
contact <- parameter()
dim(contact) <- c(n_group, n_group)

# [group, vaccination status] parameters
sus          <- parameter()   # relative susceptibility
trans        <- parameter()   # relative infectiousness when infected
p_hosp       <- parameter()   # P(hospitalised | infected)
p_death_comm <- parameter()   # P(community death | infected, not hospitalised)
p_death_hosp <- parameter()   # P(death | hospitalised)
dur_hosp     <- parameter()   # mean length of hospital stay (days)
dim(sus, trans, p_hosp, p_death_comm, p_death_hosp, dur_hosp) <- c(n_group, 2)

S_ini <- parameter()
E_ini <- parameter()
I_ini <- parameter()
R_ini <- parameter()
dim(S_ini, E_ini, I_ini, R_ini) <- c(n_group, 2)

# Piecewise-constant contact multiplier per group. Contacts between groups i
# and j are scaled by sqrt(beta_mult[i] * beta_mult[j]), so a multiplier m
# applied to everyone scales all transmission by m.
beta_time       <- parameter()
beta_mult_value <- parameter()
dim(beta_time)       <- n_beta_time
dim(beta_mult_value) <- c(n_group, n_beta_time)
beta_mult <- interpolate(beta_time, beta_mult_value, "constant")
dim(beta_mult) <- n_group

## ---- Force of infection ---------------------------------------------------
N[] <- sum(S[i, ]) + sum(E[i, ]) + sum(I[i, ]) + sum(H[i, ]) + sum(R[i, ])
dim(N) <- n_group

sqrt_mult[] <- sqrt(beta_mult[i])
dim(sqrt_mult) <- n_group

# Infectious pressure exerted by group j (per member of group j)
prev[] <- sqrt_mult[i] * (trans[i, 1] * I[i, 1] + trans[i, 2] * I[i, 2]) / N[i]
dim(prev) <- n_group

contact_prev[, ] <- contact[i, j] * prev[j]
dim(contact_prev) <- c(n_group, n_group)

lambda[] <- beta * sqrt_mult[i] * sum(contact_prev[i, ])
dim(lambda) <- n_group

## ---- Flows (binomial draws over one step of length dt) ---------------------
p_EI <- 1 - exp(-sigma * dt)
p_I_out <- 1 - exp(-gamma * dt)
p_SE[, ] <- 1 - exp(-lambda[i] * sus[i, j] * dt)
p_H_out[, ] <- 1 - exp(-dt / dur_hosp[i, j])
dim(p_SE, p_H_out) <- c(n_group, 2)

n_SE[, ] <- Binomial(S[i, j], p_SE[i, j])
n_EI[, ] <- Binomial(E[i, j], p_EI)
n_I_out[, ] <- Binomial(I[i, j], p_I_out)
n_IH[, ] <- Binomial(n_I_out[i, j], p_hosp[i, j])
n_ID[, ] <- Binomial(n_I_out[i, j] - n_IH[i, j], p_death_comm[i, j])
n_H_out[, ] <- Binomial(H[i, j], p_H_out[i, j])
n_HD[, ] <- Binomial(n_H_out[i, j], p_death_hosp[i, j])
dim(n_SE, n_EI, n_I_out, n_IH, n_ID, n_H_out, n_HD) <- c(n_group, 2)

## ---- States ---------------------------------------------------------------
initial(S[, ]) <- S_ini[i, j]
initial(E[, ]) <- E_ini[i, j]
initial(I[, ]) <- I_ini[i, j]
initial(H[, ]) <- 0
initial(R[, ]) <- R_ini[i, j]
initial(D_hosp[, ]) <- 0
initial(D_comm[, ]) <- 0
initial(cum_infections[, ]) <- 0
initial(cum_hosp[, ]) <- 0
# Daily incidence: flow during (t - 1, t]
initial(inc_infections[, ], zero_every = 1) <- 0
initial(inc_hosp[, ], zero_every = 1) <- 0
initial(inc_deaths_hosp[, ], zero_every = 1) <- 0
initial(inc_deaths_comm[, ], zero_every = 1) <- 0
dim(S, E, I, H, R, D_hosp, D_comm, cum_infections, cum_hosp,
    inc_infections, inc_hosp, inc_deaths_hosp, inc_deaths_comm) <- c(n_group, 2)

update(S[, ]) <- S[i, j] - n_SE[i, j]
update(E[, ]) <- E[i, j] + n_SE[i, j] - n_EI[i, j]
update(I[, ]) <- I[i, j] + n_EI[i, j] - n_I_out[i, j]
update(H[, ]) <- H[i, j] + n_IH[i, j] - n_H_out[i, j]
update(R[, ]) <- R[i, j] + n_I_out[i, j] - n_IH[i, j] - n_ID[i, j] +
  n_H_out[i, j] - n_HD[i, j]
update(D_hosp[, ]) <- D_hosp[i, j] + n_HD[i, j]
update(D_comm[, ]) <- D_comm[i, j] + n_ID[i, j]
update(cum_infections[, ]) <- cum_infections[i, j] + n_SE[i, j]
update(cum_hosp[, ]) <- cum_hosp[i, j] + n_IH[i, j]
update(inc_infections[, ]) <- inc_infections[i, j] + n_SE[i, j]
update(inc_hosp[, ]) <- inc_hosp[i, j] + n_IH[i, j]
update(inc_deaths_hosp[, ]) <- inc_deaths_hosp[i, j] + n_HD[i, j]
update(inc_deaths_comm[, ]) <- inc_deaths_comm[i, j] + n_ID[i, j]
