## Example 8: vaccination strategies
##
## People can be vaccinated before the epidemic (`vacc_coverage`: a fraction
## of each group starts in the vaccinated stratum) or during it
## (`vacc_schedule`: doses per day from given times, for chosen age groups or
## regions). Either way, vaccinated people keep their status for the rest of
## the run.
##
## This example looks at what each kind of vaccine effectiveness (VE) does,
## how much coverage it takes to bring R_eff below one, whether a fixed
## number of doses is better spent protecting the old directly or slowing
## transmission among the young, and how much it matters when a campaign
## starts.
##
## Run line by line (e.g. Ctrl+Enter in RStudio), or with
##   source("examples/08_vaccination.R", echo = TRUE)
## so that tables and plots are shown.

library(seirvax)
library(data.table)
library(ggplot2)

ex <- example_data()
N_age <- rowSums(ex$N)          # national population by age
N_age

base <- build_model(
  N = N_age, C_age = ex$C_age, R0 = 1.5,
  susceptibility = ex$epi$susceptibility,
  p_hosp = ex$epi$p_hosp,
  p_death_hosp = ex$epi$p_death_hosp,
  p_death_comm = ex$epi$p_death_comm,
  I0 = c(5, 10, 5, 0))
base                            # nobody vaccinated: R_eff(0) = R0

# Collect the outcomes of several deterministic scenarios in one wide table
outcomes <- function(models, t_max = 365) {
  cmp <- compare_scenarios(models, t_max = t_max)
  tab <- dcast(summary(cmp), scenario ~ measure, value.var = "value")
  tab[, R_eff := vapply(models, R_eff, 1)]
  tab[, .(scenario, R_eff = round(R_eff, 2), cum_infections = round(cum_infections),
          cum_hosp = round(cum_hosp), D = round(D))]
}

# ---- What each VE does ------------------------------------------------------
# 80% of the over-65s vaccinated, with a vaccine that has VE = 40% against one
# thing at a time. ve_hosp and ve_death are protection *given infection*, so
# they act on top of ve_infection.
elderly <- c(0, 0, 0, 0.8)
outcomes(list(
  "no vaccine"          = base,
  "VE infection 40%"    = update_model(base, vacc_coverage = elderly, ve_infection = 0.4),
  "VE transmission 40%" = update_model(base, vacc_coverage = elderly, ve_transmission = 0.4),
  "VE hosp 40%"         = update_model(base, vacc_coverage = elderly, ve_hosp = 0.4),
  "VE death 40%"        = update_model(base, vacc_coverage = elderly, ve_death = 0.4)))
# Protection against infection and transmission lowers R_eff and so helps the
# unvaccinated too; protection against hospitalisation and death only helps
# the vaccinated, and barely changes the epidemic itself. Vaccinating 65+ alone
# hardly moves R_eff, since they have few contacts, but it averts a large share
# of admissions and deaths because that is where the severe cases are.

# ---- Coverage needed to stop growth -----------------------------------------
# R_eff at the start, vaccinating the same share of every age group. With VE
# against infection e and coverage c, a homogeneous population needs
# c * e > 1 - 1/R0; the age structure shifts that a little.
grid <- CJ(coverage = seq(0, 0.95, by = 0.05), ve_infection = c(0.4, 0.6, 0.8))
grid[, R_eff := mapply(function(cv, ve) {
  R_eff(update_model(base, vacc_coverage = cv, ve_infection = ve))
}, coverage, ve_infection)]

ggplot(grid, aes(coverage, R_eff, colour = factor(ve_infection))) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x = "Coverage (all ages)", y = "R_eff at the start",
       colour = "VE against\ninfection") +
  theme_bw()

grid[R_eff < 1, .(coverage_needed = min(coverage)), by = ve_infection]

# ---- The same doses, given to different groups -------------------------------
# Enough doses for 80% of the over-65s. Give them to the over-65s, who are the
# most likely to be hospitalised, or to the 0-19s, who have the most contacts,
# or spread them over everyone. VE 40% against infection and 50% against
# hospitalisation and death given infection.
doses <- 0.8 * N_age[["65+"]]
coverage_for <- function(groups) {
  cv <- setNames(rep(0, length(N_age)), names(N_age))
  cv[groups] <- doses / sum(N_age[groups])
  cv
}
vaccine <- list(ve_infection = 0.4, ve_hosp = 0.5, ve_death = 0.5)
strategy <- function(groups) {
  do.call(update_model, c(list(base, vacc_coverage = coverage_for(groups)), vaccine))
}
strategies <- list(
  "no vaccine" = base,
  "65+"        = strategy("65+"),
  "0-19"       = strategy("0-19"),
  "everyone"   = strategy(names(N_age)))
sapply(strategies[-1], function(m) round(m$epi$vacc_coverage[, 1], 2))

strat <- outcomes(strategies)
strat
# Vaccinating the 0-19s prevents the most infections, but vaccinating the
# 65+ prevents far more admissions and deaths. Per 1,000 doses:
strat[, .(scenario,
          hosp_averted_per_1000 = round(1000 * (cum_hosp[1] - cum_hosp) / doses, 1),
          deaths_averted_per_1000 = round(1000 * (D[1] - D) / doses, 2))][-1]
# Which strategy wins depends on R0, on how steeply severity rises with age,
# and on the kind of protection: try R0 = 1.2 or 2, or ve_infection = 0, with
# update_model(base, ...) above.

# ---- Who is protected ---------------------------------------------------------
# Admissions by age and vaccination status. Vaccinating the 0-19s also
# protects the over-65s, but only indirectly, through fewer infections, and
# much less than vaccinating them directly.
cmp <- compare_scenarios(strategies[c("no vaccine", "0-19", "65+")], t_max = 365)
adm <- get_output(cmp, "cum_hosp", by = c("age", "vacc_status"))[time == 365]
dcast(adm, scenario + age ~ vacc_status, value.var = "value",
      fun.aggregate = function(x) round(sum(x)))

plot(cmp, var = "inc_hosp", by = "age")

# ---- Vaccinating during the epidemic ------------------------------------------
# vacc_schedule gives doses per day from given times. The same doses as above,
# given to the over-65s over 30 days, with the campaign starting at different
# times; "before" is the initial coverage used above, for comparison.
campaign <- function(start) {
  data.frame(time = c(start, start + 30), doses = c(doses / 30, 0), age = "65+")
}
starts <- c(0, 60, 120, 150, 180)
timing <- c(list("no vaccine" = base, "before" = strategies[["65+"]]),
            lapply(setNames(starts, paste("from day", starts)), function(t) {
              do.call(update_model, c(list(base, vacc_schedule = campaign(t)), vaccine))
            }))
timing[["from day 120"]]              # print() shows the schedule
# (R_eff is left out: it is computed at day 0, so it only sees vacc_coverage)
outcomes(timing)[, !"R_eff"]
# Without a vaccine, admissions peak around day 175. A campaign that is over
# well before the peak does almost as well as vaccinating beforehand; one
# that starts at the peak averts less than a quarter as many admissions,
# because by then most of the epidemic among the over-65s has happened.

cmp_timing <- compare_scenarios(timing, t_max = 365)
plot(cmp_timing, var = "inc_hosp")
plot(cmp_timing, var = "cum_doses")

