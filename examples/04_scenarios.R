## Example 4: comparing scenarios
##
## A typical exercise workflow: define a baseline, vary one thing at a time
## with update_model(), run each scenario stochastically and compare key
## outcomes in a table and a plot.

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

runs <- lapply(scenarios, run_model, type = "stochastic", t_max = 365,
               n_sim = 100, seed = 1, n_threads = 4)

# Key outcomes with 95% intervals
tab <- rbindlist(lapply(names(runs), function(s) summary(runs[[s]])[, scenario := s]))
dcast(tab, scenario ~ measure, value.var = "median")
tab[measure == "peak_hosp_occupancy"]

# Hospital occupancy over time
h <- rbindlist(lapply(names(runs), function(s) {
  get_output(runs[[s]], "H")[, scenario := s]
}))
h <- h[, .(median = median(value),
           lower = quantile(value, 0.025), upper = quantile(value, 0.975)),
       by = .(scenario, time)]
print(ggplot(h, aes(time, median, colour = scenario, fill = scenario)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
        geom_line() +
        geom_hline(yintercept = 400, linetype = 2) +   # e.g. hospital capacity
        labs(x = "Day", y = "Hospital occupancy") + theme_bw())
