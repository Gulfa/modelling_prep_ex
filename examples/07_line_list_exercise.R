## Example 7: generating exercise data and fitting it
##
## simulate_line_list() runs an epidemic in one of the bundled countries and
## returns a line list of cases: every hospitalisation and every death, plus a
## fraction of milder infections. This is the data an exercise would hand out.
## The second half plays the participant: build a model, fit R0 to the
## admissions, and project forward.
##
## Run line by line (e.g. Ctrl+Enter in RStudio), or with
##   source("examples/07_line_list_exercise.R", echo = TRUE)
## so that tables and plots are shown.

library(seirvax)
library(data.table)
library(ggplot2)

# ---- Organiser: make the data -------------------------------------------------
# R0 = NULL draws one at random, so the organiser does not pick it either.
# The disease is described by arguments: the defaults are a fast respiratory
# infection in an unvaccinated population. Here we slow it down a little and
# vaccinate the elderly, to show how that is done.
ll <- simulate_line_list(
  country = "Sweway",
  R0 = 1.8,                      # or NULL for a random one
  t_max = 150,
  p_detect = 0.1,                # 10% of non-severe infections are found
  age_breaks = c(0, 20, 45, 65), # how ages are reported
  latent_period = 2,
  infectious_period = 4,
  vacc_coverage = function(age) ifelse(age >= 65, 0.8, 0),
  seed = 1)

head(ll)
nrow(ll)
ll[, .(cases = .N, hospitalised = sum(hospitalised), deaths = sum(died)),
   by = age_group]

# The truth is attached, for checking answers afterwards:
truth <- attr(ll, "truth")
truth$R0
truth$model                      # the model that generated the data

# ---- Participant: what the data looks like ------------------------------------
# Epidemic curve by age group
curve <- ll[, .N, by = .(date_onset, age_group)]
ggplot(curve, aes(date_onset, N, fill = age_group)) +
  geom_col() +
  labs(x = NULL, y = "Cases by onset date") + theme_bw()

# Admissions and deaths over time
severe <- rbind(
  ll[!is.na(date_hospital), .(date = date_hospital, event = "admission")],
  ll[!is.na(date_death), .(date = date_death, event = "death")])
ggplot(severe[, .N, by = .(date, event)], aes(date, N, colour = event)) +
  geom_line() + labs(x = NULL, y = "Per day") + theme_bw()

# Which region started it, and how severe are the cases by age?
ll[, .(first_case = min(date_onset), cases = .N), by = region][order(first_case)]
ll[, .(hosp_share = mean(hospitalised), cfr = mean(died)), by = age_group]

# ---- Participant: fit R0 ------------------------------------------------------
# Daily admissions in the format fit_model() wants. Day 1 of the model is the
# first onset in the data, so the epidemic starts when the data does.
start <- min(ll$date_onset)
data <- admissions_from_line_list(ll, start_date = start)
data <- data[time <= 60]         # only the first 60 days are "available" yet
head(data)

# A model of the same country, at the age resolution we report. The participant
# does not know R0, the seeding or the true severity.
N <- population_matrix("Sweway", age_breaks = c(0, 20, 45, 65))

# Regions mix by distance (see example 06 and docs/data.md)
dist <- region_distances(country = "Sweway")
w <- exp(-dist / 100); diag(w) <- 0
M <- 0.1 * w / rowSums(w); diag(M) <- 0.9

# Seed where the first cases were seen
first_region <- ll[, .(first = min(date_onset)), by = region][order(first)][1]$region
I0 <- matrix(0, nrow(N), ncol(N), dimnames = dimnames(N))
I0["20-44", first_region] <- 20

# The natural history and severity are what a participant would assume from
# the literature (and from the line list itself) -- not necessarily what
# generated the data.
model <- build_model(
  N = N,
  C_age = age_contact_matrix(rownames(N)),
  M_region = M,
  R0 = 2.5,                      # a starting guess; this is what we fit
  latent_period = 2,
  infectious_period = 4,
  susceptibility = c(0.7, 1, 1, 1),
  p_hosp = c(0.01, 0.03, 0.1, 0.3),
  p_death_hosp = c(0.1, 0.15, 0.25, 0.4),
  vacc_coverage = c(0, 0, 0, 0.8),
  ve_infection = 0.4, ve_hosp = 0.7, ve_death = 0.5,
  I0 = I0)

# The age pattern of severity in the data is a useful check on those guesses:
ll[, .(hospitalised_share = round(mean(hospitalised), 3)), by = age_group]

fit <- fit_model(model, data, pars = "R0", method = "mcmc",
                 n_steps = 2000, n_chains = 4, seed = 1)
fit
truth$R0                         # how close was it?
# Expect it to be close but not exact, and the credible interval to be too
# narrow to cover the truth: the fitted model has four age groups rather than
# twenty, assumes its own severity, and treats one stochastic epidemic as if
# it were deterministic. The interval measures uncertainty *within* the model,
# not the error of the model itself.
plot(fit, t_max = 150)
plot(fit, type = "trace")

# Fitting a severity multiplier as well absorbs the fact that we do not know
# the seeding or the true severity; with only the early epidemic the two are
# hard to separate (see docs/model_description.md).
fit2 <- fit_model(model, data, pars = c("R0", "hosp_scale"))
fit2

# ---- Participant: project and compare scenarios --------------------------------
# project_fit() starts where the data ends, so this is the forecast proper.
# For a pmcmc fit each draw continues from the particle filter's state, so it
# carries on the epidemic that was observed rather than an average one.
proj <- project_fit(fit, t_max = 150, n_draws = 100, vars = "H")
proj[, .(median = median(value), upper = quantile(value, 0.975)), by = time][
  time %in% c(60, 90, 120, 150)]

# Who ends up in hospital? Totals by age and region with the fit's uncertainty,
# and each stratum's share of the national total.
project_summary(fit, by = "age", t_max = 365, n_draws = 100, seed = 1)
project_summary(fit, by = "region", t_max = 365, n_draws = 100,
                measures = c("cum_hosp", "D"), seed = 1)
plot(fit, type = "distribution", by = "age", measure = "cum_hosp",
     t_max = 365, n_draws = 100, seed = 1)

# Note how narrow the age shares are when only R0 is fitted: the age pattern
# comes from the contact matrix and the severity assumptions, not from the data.

# What actually happened, from the full data:
admissions_from_line_list(ll, start_date = start)[time %in% c(60, 90, 120, 150)]

# A scenario on top of the fit: 30% contact reduction from day 60. Passing the
# scenario model to project_fit()/project_summary() keeps the uncertainty from
# the fit; update_model(fit$model, ...) would give a single trajectory at the
# point estimate instead.
lockdown <- update_model(fit$base_model,
                         beta_schedule = data.frame(time = 60, mult = 0.7))
project_summary(fit, t_max = 200, n_draws = 100, measures = "cum_hosp", seed = 1)
project_summary(fit, model = lockdown, t_max = 200, n_draws = 100,
                measures = "cum_hosp", seed = 1)
