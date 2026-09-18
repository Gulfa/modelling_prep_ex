## Example 1: a single, homogeneous population
##
## The simplest use of the model: one group, no age or region structure.
## Shows how to build a model, run the deterministic and stochastic versions,
## and look at the results.

library(seirvax)

# ---- Build the model --------------------------------------------------------
# All parameters are set here; nothing is simulated yet.
m <- build_model(
  N = 500000,               # population
  R0 = 1.8,                 # beta is calibrated from R0
  latent_period = 3,        # days
  infectious_period = 5,    # days
  p_hosp = 0.02,            # probability an infection is hospitalised
  p_death_hosp = 0.1,       # probability a hospitalised case dies
  p_death_comm = 0.0005,    # probability a non-hospitalised case dies
  dur_hosp = 7,             # mean length of stay (days)
  I0 = 20)                  # initial infectious people
m

# ---- Deterministic run ------------------------------------------------------
det <- run_model(m, type = "deterministic", t_max = 250)
summary(det)
print(plot(det, var = "inc_infections"))
print(plot(det, var = "H"))           # hospital occupancy

# ---- Stochastic run ---------------------------------------------------------
sto <- run_model(m, type = "stochastic", t_max = 250, n_sim = 100, seed = 1)
summary(sto)                          # mean, median and 95% interval
print(plot(sto, var = "H", n_traces = 20))

# ---- Working with the output directly ----------------------------------------
# get_output() returns a data.table in long format: time, sim, variable, value
d <- get_output(sto, vars = c("inc_hosp", "H"))
head(d)

# e.g. the probability that hospital occupancy ever exceeds 300 beds
peak <- d[variable == "H", .(peak = max(value)), by = sim]
mean(peak$peak > 300)
