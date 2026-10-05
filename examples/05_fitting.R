## Example 5: fitting R0 to daily hospital admissions
##
## We simulate "observed" admissions from a model with a known R0, then fit
## R0 back from the first 60 days, and project forward. The three fitting
## methods are shown: a quick point estimate ("optim"), Bayesian MCMC on the
## deterministic model ("mcmc"), and particle MCMC on the stochastic model
## ("pmcmc").
##
## Run line by line (e.g. Ctrl+Enter in RStudio), or with
##   source("examples/05_fitting.R", echo = TRUE)
## so that tables and plots are shown.

library(seirvax)
library(data.table)

ex <- example_data()
I0 <- matrix(0, 4, 3, dimnames = dimnames(ex$N))
I0["20-44", ] <- 10

# ---- A "true" epidemic and some data ----------------------------------------
truth <- build_model(
  N = ex$N, C_age = ex$C_age, M_region = ex$M_region,
  R0 = 1.8,
  susceptibility = ex$epi$susceptibility,
  p_hosp = ex$epi$p_hosp,
  p_death_hosp = ex$epi$p_death_hosp,
  vacc_coverage = ex$epi$vacc_coverage,
  ve_infection = 0.4, ve_hosp = 0.7,
  I0 = I0)

# One stochastic realisation plays the role of the observed data: daily
# hospital admissions for the first 60 days. The data needs columns `time`
# (day since the model start) and `admissions`.
sim <- run_model(truth, "stochastic", t_max = 60, n_sim = 1, seed = 42)
data <- get_output(sim, "inc_hosp")[time > 0, .(time, admissions = value)]
tail(data)

# ---- The model to fit ---------------------------------------------------------
# Everything except the fitted parameters is taken from this model. Here we
# pretend not to know R0 and start from a wrong guess.
model <- update_model(truth, R0 = 2.5)
fit_parameters()          # what can be fitted, with default priors

# ---- 1. Point estimate (fast) -----------------------------------------------
fit1 <- fit_model(model, data, pars = "R0", method = "optim")
fit1
plot(fit1, t_max = 150)   # data, fit and projection with uncertainty

# ---- 2. Bayesian MCMC on the deterministic model ------------------------------
# Also estimates the dispersion of the data (nb_size). Priors can be changed:
fit2 <- fit_model(model, data, pars = c("R0", "nb_size"), method = "mcmc",
                  prior = c(R0 = "Uniform(1, 4)"),
                  n_steps = 2000, n_chains = 4, seed = 1)
fit2
plot(fit2, type = "trace")   # chains should overlap and look like noise
plot(fit2, t_max = 150)

# ---- 3. Particle MCMC on the stochastic model ---------------------------------
# Accounts for randomness in the epidemic itself, but each step runs a
# particle filter, so it is much slower: keep the chain short here.
fit3 <- fit_model(model, data, pars = "R0", method = "pmcmc",
                  n_steps = 300, n_chains = 1, n_particles = 100, dt = 0.25,
                  n_threads = 4, seed = 1)
fit3

# ---- Compare the estimates -----------------------------------------------------
rbind(cbind(method = "optim", summary(fit1)),
      cbind(method = "mcmc", summary(fit2)),
      cbind(method = "pmcmc", summary(fit3)))

# ---- Using the fit -------------------------------------------------------------
# fit$model is the model at the estimate and works like any other model:
summary(run_model(fit2$model, t_max = 300))

# Projections over posterior draws, e.g. hospital occupancy by age group:
# project_fit() returns a run-like object covering the days after the data,
# with one simulation per posterior draw.
proj <- project_fit(fit2, t_max = 200, n_draws = 50, seed = 1)
proj
plot(proj, var = "H", by = "age")
get_output(proj, "H", by = "age")[time == 120,
  .(median = median(value), upper = quantile(value, 0.975)), by = age]

# Scenarios can be built on the fitted model with update_model(), e.g. a
# 30% contact reduction from day 70:
lockdown <- update_model(fit2$model,
                         beta_schedule = data.frame(time = 70, mult = 0.7))
plot(run_model(lockdown, t_max = 250), var = "inc_hosp")

# ---- Admissions by region ------------------------------------------------------
# With `region` (and/or `age`) columns, each region is its own data series.
data_region <- get_output(sim, "inc_hosp", by = "region")[
  time > 0, .(time, region, admissions = value)]
fit4 <- fit_model(model, data_region, pars = "R0")
fit4
plot(fit4, t_max = 120)
