# seirvax

An SEIR model with hospitalisation, death and vaccination, built for preparedness exercises.

- Stratified by **age group** and **region**. Age and region contact matrices are combined into one contact structure.
- Has a **vaccinated stratum** that runs in parallel with the unvaccinated one. Vaccine effectiveness against infection, onward transmission, hospitalisation and death can vary by age.
- Hospitalisation, **hospital deaths** and **community deaths**.
- **Contact multipliers over time** (interventions), either global or targeted at specific regions or age groups, plus **mixing matrices that can change over time** (school holidays, travel restrictions).
- A **deterministic** (ODE) and a **stochastic** (tau-leaping) version of the same model, written in [odin2](https://mrc-ide.github.io/odin2/) and run with [dust2](https://mrc-ide.github.io/dust2/).
- **Bundled data** for four fictional countries: population by region and five-year age band, plus region polygons, settlements and a population grid for maps.
- **Synthetic exercise data**: a line list of cases with hospitalisation and death dates, generated from a real run of the model.
- **Fitting to daily hospital admissions** (national, or by region or age): a quick point estimate, MCMC, or particle MCMC, using [monty](https://mrc-ide.github.io/monty/).

## Installation

The package needs `dust2` and `monty` from the mrc-ide r-universe, and a C++ compiler to build it.

```r
install.packages(c("dust2", "monty", "data.table", "ggplot2", "remotes"),
                 repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
remotes::install_local("path/to/modelling_prep_ex")
```

For working on the model itself, see [`docs/modifying_the_model.md`](docs/modifying_the_model.md#development-setup).

## Quick start

```r
library(seirvax)

ex <- example_data()                      # made-up example inputs (4 ages x 3 regions)

m <- build_model(
  N = ex$N, C_age = ex$C_age, M_region = ex$M_region,
  R0 = 2,
  susceptibility = ex$epi$susceptibility,
  p_hosp = ex$epi$p_hosp,
  p_death_hosp = ex$epi$p_death_hosp,
  p_death_comm = ex$epi$p_death_comm,
  vacc_coverage = ex$epi$vacc_coverage,
  ve_infection = 0.4, ve_hosp = 0.7, ve_death = 0.5,
  I0 = 5,
  beta_schedule = data.frame(time = c(60, 120), mult = c(0.6, 1)))
m                                          # parameters, R0 and R_eff(0)

det <- run_model(m, "deterministic", t_max = 300)
sto <- run_model(m, "stochastic", t_max = 300, n_sim = 100, seed = 1)

summary(sto, by = "age")                   # cum_infections, cum_hosp, D, peak_hosp_*
plot(sto, var = "H", by = "region")        # hospital occupancy by region
get_output(det, vars = "inc_hosp", by = c("age", "vacc_status"))   # tidy data.table

m2 <- update_model(m, ve_hosp = 0.5)       # a scenario variant
```

Fitting R0 to daily admissions (a data.frame with columns `time` and `admissions`):

```r
fit <- fit_model(m, data, pars = "R0", method = "mcmc")
fit                                        # estimates with 95% intervals
plot(fit, t_max = 200)                     # data, fit and projection
run_model(fit$model, t_max = 300)          # the fitted model, for scenarios
```

## Main functions

| Function | Purpose |
|---|---|
| `build_model()` | Set all parameters; returns a model object (nothing is simulated) |
| `update_model()` | Copy of a model with some parameters changed |
| `run_model()` | Run the `"deterministic"` or `"stochastic"` version |
| `get_output()` | Results as a long `data.table`, summed to the strata you choose (`by`) |
| `summary()` | Totals and peaks (with intervals for stochastic runs), named as in `get_output()` |
| `plot()` | Plot one variable over time; stochastic runs as a spaghetti plot, or `style = "ribbon"` |
| `R0()`, `R_eff()` | Reproduction numbers from the next-generation matrix |
| `fit_model()` | Fit `R0` (and optionally `hosp_scale`, `nb_size`) to daily hospital admissions |
| `fit_parameters()` | Parameters that can be fitted, with default priors |
| `project_fit()` | Run the model over posterior draws, e.g. to project beyond the data |
| `project_summary()` | Totals and peaks by age and region from a fit, with uncertainty |
| `compare_scenarios()` | Several models, runs or projections side by side, with `summary()` and `plot()` |
| `example_data()` | Made-up example population, contact matrices and parameters |
| `read_demography()`, `region_info()` | The bundled country data, and one row per region |
| `population_matrix()` | Age × region population matrix for `build_model(N = ...)` |
| `read_geography()`, `region_distances()`, `plot_region_map()` | Map layers, distances between region centres, and maps of model output |
| `simulate_line_list()` | Generate a line list of cases for an exercise, with the truth attached |
| `admissions_from_line_list()` | Daily admissions from a line list, ready for `fit_model()` |
| `age_contact_matrix()` | A synthetic contact matrix for any set of age groups |

Parameters can be given as a single value, a vector by age group, a vector by region, or an age × region matrix. See `?build_model`.

## Examples

The scripts in [`examples/`](examples/) are the best place to start:

1. [`01_single_population.R`](examples/01_single_population.R): one group, deterministic vs. stochastic.
2. [`02_age_structured.R`](examples/02_age_structured.R): age mixing, susceptibility, severity and vaccine effectiveness by age.
3. [`03_age_and_region.R`](examples/03_age_and_region.R): combined age and region mixing, regional seeding, global and targeted interventions.
4. [`04_scenarios.R`](examples/04_scenarios.R): comparing several scenarios with `compare_scenarios()`.
5. [`05_fitting.R`](examples/05_fitting.R): fitting R0 to hospital admissions with each method, and projecting from the fit.
6. [`06_country_data.R`](examples/06_country_data.R): building a model for one of the bundled countries, and mapping the results.
7. [`07_line_list_exercise.R`](examples/07_line_list_exercise.R): generating a line list for an exercise, then fitting R0 to it and projecting.

## Documentation

- [`docs/data.md`](docs/data.md): the bundled country data and how to turn it into model inputs.
- [`docs/model_description.md`](docs/model_description.md): model structure, equations, assumptions and known limitations.
- [`docs/modifying_the_model.md`](docs/modifying_the_model.md): how the code is organised, and how to change the model, with a worked example.
- Function help: `?build_model`, `?run_model`, `?get_output`, and so on.
