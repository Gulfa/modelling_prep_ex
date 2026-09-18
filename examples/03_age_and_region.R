## Example 3: age groups and regions
##
## The population is split by age and region. The age contact matrix C_age is
## combined with a region mixing matrix M_region, where M_region[r, s] is the
## fraction of a region-r resident's contacts made with residents of region s.
## The epidemic is seeded in one region and a regional lockdown is added.

library(seirvax)
library(data.table)

ex <- example_data()
ex$N           # age x region population
ex$M_region    # rows sum to 1

# Seed 20 infections among 20-44 year olds in the North
I0 <- matrix(0, 4, 3, dimnames = dimnames(ex$N))
I0["20-44", "North"] <- 20

m <- build_model(
  N = ex$N,
  C_age = ex$C_age,
  M_region = ex$M_region,
  R0 = 2,
  susceptibility = ex$epi$susceptibility,
  p_hosp = ex$epi$p_hosp,
  p_death_hosp = ex$epi$p_death_hosp,
  p_death_comm = ex$epi$p_death_comm,
  # Parameters can also be age x region matrices, e.g. lower coverage in South
  vacc_coverage = cbind(North = ex$epi$vacc_coverage,
                        Central = ex$epi$vacc_coverage,
                        South = ex$epi$vacc_coverage * 0.7),
  ve_infection = 0.4, ve_hosp = 0.7, ve_death = 0.5,
  I0 = I0)
m

res <- run_model(m, t_max = 300)
print(plot(res, var = "inc_infections", by = "region"))
summary(res, by = "region")

# ---- Contact multipliers (beta_schedule) -------------------------------------
# A global schedule: contacts reduced by 40% for everyone from day 60 to 120.
global <- data.frame(time = c(60, 120), mult = c(0.6, 1))

# Targeted: rows with `region` and/or `age` only apply to those groups. Rows
# are applied in time order and override earlier ones for the matching groups.
# Here: the global measure, plus a stricter lockdown in North from day 40
# to day 120, and older people reducing contacts for the whole period.
targeted <- data.frame(
  time   = c(40,      60,  60,    120,  120),
  mult   = c(0.4,     0.6, 0.5,   1,    0.8),
  region = c("North", NA,  NA,    NA,   NA),
  age    = c(NA,      NA,  "65+", NA,   "65+"))
targeted
# Note: at day 60 the global row also applies to North, ending its stricter
# lockdown. To keep North at 0.4, add a row (time = 60, mult = 0.4,
# region = "North") after the global one.

res_global <- run_model(update_model(m, beta_schedule = global), t_max = 300)
res_target <- run_model(update_model(m, beta_schedule = targeted), t_max = 300)

h <- rbind(
  get_output(res, "H", by = "region")[, scenario := "no measures"],
  get_output(res_global, "H", by = "region")[, scenario := "global"],
  get_output(res_target, "H", by = "region")[, scenario := "targeted"])
library(ggplot2)
print(ggplot(h, aes(time, value, colour = scenario)) +
        geom_line() + facet_wrap(~region, scales = "free_y") +
        labs(x = "Day", y = "Hospital occupancy") + theme_bw())

# ---- Stochastic version: early fade-out is possible --------------------------
sto <- run_model(m, "stochastic", t_max = 300, n_sim = 100, seed = 1, n_threads = 4)
print(plot(sto, var = "H", by = "region", n_traces = 20))
final <- get_output(sto, "cum_infections")[time == 300]
mean(final$value < 1000)      # fraction of simulations where the epidemic died out
