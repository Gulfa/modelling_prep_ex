## Example 4: comparing scenarios
##
## A typical exercise workflow: define a baseline, vary one thing at a time
## with update_model(), run each scenario stochastically and compare key
## outcomes in a table and a plot.
##
## Run line by line (e.g. Ctrl+Enter in RStudio), or with
##   source("examples/04_scenarios.R", echo = TRUE)
## so that tables and plots are shown.

library(seirvax)
library(data.table)
library(ggplot2)

ex <- example_data()
I0 <- matrix(0, 4, 3, dimnames = dimnames(ex$N))
I0["20-44", ] <- 10

baseline <- build_model(
  N = ex$N, C_age = ex$C_age, M_region = ex$M_region,
  R0 = 1.8,
  susceptibility = ex$epi$susceptibility,
  p_hosp = ex$epi$p_hosp,
  p_death_hosp = ex$epi$p_death_hosp,
  p_death_comm = ex$epi$p_death_comm,
  vacc_coverage = ex$epi$vacc_coverage,
  ve_infection = 0.4, ve_hosp = 0.7, ve_death = 0.5,
  I0 = I0)

scenarios <- list(
  "baseline"          = baseline,
  "no vaccine"        = update_model(baseline, vacc_coverage = 0),
  "weaker vaccine"    = update_model(baseline, ve_infection = 0.2, ve_hosp = 0.5),
  "higher R0"         = update_model(baseline, R0 = 2.3),
  "distancing day 50" = update_model(baseline,
                          beta_schedule = data.frame(time = c(50, 110), mult = c(0.7, 1))))

# compare_scenarios() runs them all and keeps them together. Models are run
# for you; runs from run_model() and projections from project_fit() can be
# passed in directly instead.
cmp <- compare_scenarios(scenarios, type = "stochastic", t_max = 365,
                         n_sim = 100, seed = 1, n_threads = 4)
cmp

# Key outcomes with 95% intervals
summary(cmp)
dcast(summary(cmp), scenario ~ measure, value.var = "median")
summary(cmp)[measure == "peak_hosp_occupancy"]
summary(cmp, by = "age")[measure == "cum_hosp"]

# Hospital occupancy over time: median and 95% interval per scenario
plot(cmp, var = "H") +
  geom_hline(yintercept = 400, linetype = 2)   # e.g. hospital capacity

plot(cmp, var = "inc_hosp", by = "region")     # one panel per region
plot(cmp, var = "H", style = "spaghetti")      # individual simulations

# get_output() gives the same long format as for a single run, with a
# scenario column, if you would rather do your own summaries
get_output(cmp, "D", by = "age")[time == 365,
  .(deaths = median(value)), by = .(scenario, age)]
