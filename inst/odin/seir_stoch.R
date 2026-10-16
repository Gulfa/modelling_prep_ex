## Stochastic SEIR-HD model with a parallel vaccinated stratum (tau-leaping).
##
## Discrete-time version of seir_ode.R: the same parameters and states, with
## binomial draws over a step of `dt` days (set when the system is created).
## Outflows from S, I and H are split with chained binomials, so the competing
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

# contact[i, j]: contacts per day a person in group i has with group j.
# Piecewise constant in time, so the mixing pattern can change (school
# holidays, travel restrictions); the last dimension of contact_value is time.
n_contact_time <- parameter(type = "integer", constant = TRUE)
contact_time   <- parameter()
contact_value  <- parameter()
dim(contact_time)  <- n_contact_time
dim(contact_value) <- c(n_group, n_group, n_contact_time)
contact <- interpolate(contact_time, contact_value, "constant")
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

# Piecewise-constant contact multiplier per group: the fraction of the group's
# contacts that still happen. Contacts between groups i and j are scaled by
# min(beta_mult[i], beta_mult[j]) -- the more restricted group decides -- so a
# multiplier m applied to everyone scales all transmission by m.
beta_time       <- parameter()
beta_mult_value <- parameter()
dim(beta_time)       <- n_beta_time
dim(beta_mult_value) <- c(n_group, n_beta_time)
beta_mult <- interpolate(beta_time, beta_mult_value, "constant")
dim(beta_mult) <- n_group

# Vaccination over time: doses per day given to each group, piecewise constant
# (built in R from vacc_schedule). Doses go to unvaccinated people who are not
# ill -- S and R, in proportion -- and move them to the vaccinated stratum.
n_vacc_time      <- parameter(type = "integer", constant = TRUE)
vacc_time        <- parameter()
vacc_doses_value <- parameter()
dim(vacc_time)        <- n_vacc_time
dim(vacc_doses_value) <- c(n_group, n_vacc_time)
vacc_doses <- interpolate(vacc_time, vacc_doses_value, "constant")
dim(vacc_doses) <- n_group

## ---- Force of infection ---------------------------------------------------
N[] <- sum(S[i, ]) + sum(E[i, ]) + sum(I[i, ]) + sum(H[i, ]) + sum(R[i, ])
dim(N) <- n_group

# Infectious pressure exerted by group j (per member of group j)
prev[] <- (trans[i, 1] * I[i, 1] + trans[i, 2] * I[i, 2]) / N[i]
dim(prev) <- n_group

contact_prev[, ] <- contact[i, j] * min(beta_mult[i], beta_mult[j]) * prev[j]
dim(contact_prev) <- c(n_group, n_group)

lambda[] <- beta * sum(contact_prev[i, ])
dim(lambda) <- n_group

## ---- Vaccination ----------------------------------------------------------
# Per-capita rate among the eligible that delivers vacc_doses per day. It is
# capped at 10 per day, so that a group whose eligible people run out is
# vaccinated to near-exhaustion within hours rather than going negative.
vacc_eligible[] <- S[i, 1] + R[i, 1]
vacc_rate[] <- min(vacc_doses[i] / max(vacc_eligible[i], 1e-10), 10)
dim(vacc_eligible, vacc_rate) <- n_group
# by stratum: only the unvaccinated are vaccinated
vacc_rate_s[, 1] <- vacc_rate[i]
vacc_rate_s[, 2] <- 0
dim(vacc_rate_s) <- c(n_group, 2)

## ---- Flows (binomial draws over one step of length dt) ---------------------
p_EI <- 1 - exp(-sigma * dt)
p_I_out <- 1 - exp(-gamma * dt)
# S is left by infection or vaccination, split with a chained binomial
rate_SE[, ] <- lambda[i] * sus[i, j]
p_S_out[, ] <- 1 - exp(-(rate_SE[i, j] + vacc_rate_s[i, j]) * dt)
p_SE_given_out[, ] <- rate_SE[i, j] / max(rate_SE[i, j] + vacc_rate_s[i, j], 1e-30)
p_RV[, ] <- 1 - exp(-vacc_rate_s[i, j] * dt)
p_H_out[, ] <- 1 - exp(-dt / dur_hosp[i, j])
dim(rate_SE, p_S_out, p_SE_given_out, p_RV, p_H_out) <- c(n_group, 2)

n_S_out[, ] <- Binomial(S[i, j], p_S_out[i, j])
n_SE[, ] <- Binomial(n_S_out[i, j], p_SE_given_out[i, j])
n_SV[, ] <- n_S_out[i, j] - n_SE[i, j]
n_RV[, ] <- Binomial(R[i, j], p_RV[i, j])
# arriving in the vaccinated stratum
n_SV_in[, 1] <- 0
n_SV_in[, 2] <- n_SV[i, 1]
n_RV_in[, 1] <- 0
n_RV_in[, 2] <- n_RV[i, 1]
dim(n_S_out, n_SV, n_RV, n_SV_in, n_RV_in) <- c(n_group, 2)
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
# Doses given, recorded in the vaccinated column (the unvaccinated one stays 0)
initial(cum_doses[, ]) <- 0
initial(inc_doses[, ], zero_every = 1) <- 0
dim(S, E, I, H, R, D_hosp, D_comm, cum_infections, cum_hosp,
    inc_infections, inc_hosp, inc_deaths_hosp, inc_deaths_comm,
    cum_doses, inc_doses) <- c(n_group, 2)

update(S[, ]) <- S[i, j] - n_S_out[i, j] + n_SV_in[i, j]
update(E[, ]) <- E[i, j] + n_SE[i, j] - n_EI[i, j]
update(I[, ]) <- I[i, j] + n_EI[i, j] - n_I_out[i, j]
update(H[, ]) <- H[i, j] + n_IH[i, j] - n_H_out[i, j]
update(R[, ]) <- R[i, j] + n_I_out[i, j] - n_IH[i, j] - n_ID[i, j] +
  n_H_out[i, j] - n_HD[i, j] - n_RV[i, j] + n_RV_in[i, j]
update(D_hosp[, ]) <- D_hosp[i, j] + n_HD[i, j]
update(D_comm[, ]) <- D_comm[i, j] + n_ID[i, j]
update(cum_infections[, ]) <- cum_infections[i, j] + n_SE[i, j]
update(cum_hosp[, ]) <- cum_hosp[i, j] + n_IH[i, j]
update(inc_infections[, ]) <- inc_infections[i, j] + n_SE[i, j]
update(inc_hosp[, ]) <- inc_hosp[i, j] + n_IH[i, j]
update(inc_deaths_hosp[, ]) <- inc_deaths_hosp[i, j] + n_HD[i, j]
update(inc_deaths_comm[, ]) <- inc_deaths_comm[i, j] + n_ID[i, j]
update(cum_doses[, ]) <- cum_doses[i, j] + n_SV_in[i, j] + n_RV_in[i, j]
update(inc_doses[, ]) <- inc_doses[i, j] + n_SV_in[i, j] + n_RV_in[i, j]

## ---- Comparison to data (used only when fitting) --------------------------
# Observed daily hospital admissions, as n_obs series. obs_map[k, g] = 1 if
# group g contributes to series k (built in R from the data's region / age
# columns; by default one national series).
n_obs  <- parameter(type = "integer", constant = TRUE)
obs_map <- parameter()
dim(obs_map) <- c(n_obs, n_group)
nb_size <- parameter()   # negative binomial dispersion (larger = less noise)

adm_contrib[, ] <- obs_map[i, j] * (inc_hosp[j, 1] + inc_hosp[j, 2])
dim(adm_contrib) <- c(n_obs, n_group)
model_admissions[] <- sum(adm_contrib[i, ])
dim(model_admissions) <- n_obs

admissions <- data()
dim(admissions) <- n_obs
admissions[] ~ NegativeBinomial(size = nb_size, mu = model_admissions[i] + 1e-6)
