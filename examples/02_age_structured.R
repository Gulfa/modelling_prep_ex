## Example 2: age structure and vaccination
##
## Four age groups with an age contact matrix, age-specific susceptibility and
## severity, and initial vaccination coverage with age-specific VE.

library(seirvax)
library(data.table)

ex <- example_data()            # made-up example inputs
N_age <- rowSums(ex$N)          # national population by age
ex$C_age                        # contacts per day by age
ex$epi                          # age-specific parameters

m <- build_model(
  N = N_age,
  C_age = ex$C_age,
  R0 = 2,
  susceptibility = ex$epi$susceptibility,   # vectors are by age group
  p_hosp = ex$epi$p_hosp,
  p_death_hosp = ex$epi$p_death_hosp,
  p_death_comm = ex$epi$p_death_comm,
  vacc_coverage = ex$epi$vacc_coverage,
  ve_infection = 0.4,                       # a single value applies to all ages
  ve_hosp = c(0.8, 0.8, 0.7, 0.6),          # ... or give one value per age group
  ve_death = 0.5,
  I0 = c(5, 10, 5, 0))
m                                           # note R_eff(0) < R0 due to vaccination

res <- run_model(m, t_max = 300)
summary(res, by = "age")
print(plot(res, var = "inc_hosp", by = "age"))
print(plot(res, var = "H", by = c("vacc_status", "age")))

# Who is hospitalised? Admissions by age and vaccination status
adm <- get_output(res, "cum_hosp", by = c("age", "vacc_status"))[time == 300]
dcast(adm, age ~ vacc_status, value.var = "value")

# Children twice as infectious (transmissibility), everything else the same
m_trans <- update_model(m, transmissibility = c(2, 1, 1, 1))
m_trans                                    # beta is recalibrated to keep R0 = 2
summary(run_model(m_trans, t_max = 300), by = "age")
