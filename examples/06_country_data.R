## Example 6: using the bundled country data
##
## The package ships fictional demographic and map data for four countries
## (Sweway, Finmark, Norland, Denden), each with five regions and five-year
## age bands. This example builds a metapopulation model for one country
## straight from that data, runs it, and maps the result.
##
## Run line by line (e.g. Ctrl+Enter in RStudio), or with
##   source("examples/06_country_data.R", echo = TRUE)
## so that tables and plots are shown.

library(seirvax)
library(data.table)

# ---- What is in the data ------------------------------------------------------
d <- read_demography()
head(d)
region_info()                      # area, density, population, settlement profile
region_info(country = "Norland")

# ---- Population for a model ---------------------------------------------------
# population_matrix() returns the age x region matrix that build_model() wants.
# age_breaks collapses the five-year bands into the groups you want to model.
N <- population_matrix("Sweway", age_breaks = c(0, 20, 45, 65))
N
colSums(N)                         # population by region
round(prop.table(N, 2), 3)         # age composition differs between regions

# Other selections:
population_matrix(regions = c("S1", "D1"), age_breaks = c(0, 65))  # across countries
population_matrix(by_region = FALSE, age_breaks = c(0, 20, 45, 65))  # whole countries

# ---- A region mixing matrix ---------------------------------------------------
# M_region[r, s] is the share of a region-r resident's contacts made in region s.
# Here we build one from the distances between region centres: most contacts at
# home, the rest spread with distance. Replace this with real travel data when
# you have it.
dist <- region_distances(country = "Sweway")
round(dist)

within <- 0.9                      # share of contacts made in the home region
w <- exp(-dist / 100)              # 100 km distance scale
diag(w) <- 0
M <- (1 - within) * w / rowSums(w)
diag(M) <- within
round(M, 3)

# ---- Build and run ------------------------------------------------------------
C_age <- matrix(c(8, 2.5, 1.5, 0.7,                 # a made-up contact matrix
                  3.0, 6, 3, 1.2,
                  1.5, 2.5, 4, 1.5,
                  0.5, 0.7, 1, 2.5), nrow = 4)

I0 <- matrix(0, nrow(N), ncol(N), dimnames = dimnames(N))
I0["20-44", "S2"] <- 20            # seeded in one region

m <- build_model(
  N = N, C_age = C_age, M_region = M,
  R0 = 1.8,
  susceptibility = c(0.7, 1, 1, 1),
  p_hosp = c(0.002, 0.01, 0.03, 0.12),
  p_death_hosp = c(0.01, 0.03, 0.08, 0.25),
  vacc_coverage = c(0, 0.3, 0.6, 0.9),
  ve_infection = 0.4, ve_hosp = 0.7, ve_death = 0.5,
  I0 = I0)
m

res <- run_model(m, t_max = 365)
summary(res, by = "region")
plot(res, var = "inc_hosp", by = "region")

# ---- Map the result -----------------------------------------------------------
attack <- get_output(res, "cum_infections", by = "region")[time == 365]
attack[, per_1000 := 1000 * value / colSums(N)[as.character(region)]]
plot_region_map(attack[, .(region, per_1000)], country = "Sweway",
                label = "Infections per 1000")

# The map layers are plain sf objects, so you can use them directly:
settlements <- read_geography("settlements", country = "Sweway")
head(settlements[order(-settlements$pop), c("name", "region", "role", "pop")])
