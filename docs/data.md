# Bundled country data

The package ships fictional demographic and map data for four countries used in the exercise. Nothing in it is real.

| | |
|---|---|
| Countries | Sweway (regions S1–S5), Finmark (F1–F5), Norland (N1–N5), Denden (D1–D5) |
| Age bands | 20 five-year bands, `0-4` … `90-94`, `95+` |
| Populations | about 0.1–3.5 million per region |
| Files | `inst/extdata/region_demographics_model_long.csv`, `inst/extdata/geography_model.gpkg` |

Region codes are unique across countries, so regions from different countries can be mixed in one model.

## Demographics

`read_demography()` returns one row per region and age band:

| Column | Meaning |
|---|---|
| `country`, `region` | country name and region code |
| `age_group` | five-year band, an **ordered** factor (so `0-4 < 5-9 < … < 95+`, not alphabetical) |
| `population` | people in that region and age band |
| `share_of_region` | that band's share of the region |
| `area_km2`, `density_per_km2`, `pop_total` | region-level values, repeated on every row |
| `profile` | settlement profile: `Metro-family`, `Urban-young`, `Mixed`, `Rural-mature`, `Remote-aging` |

`region_info()` gives the region-level columns, one row per region.

## Populations for a model

`population_matrix()` produces the age × region matrix that `build_model()` takes as `N`:

```r
N <- population_matrix("Sweway", age_breaks = c(0, 20, 45, 65))
#          S1     S2    S3     S4     S5
# 0-19  38971 134128 27538 309870 210102
# 20-44 48367 181661 28344 588507 278650
# 45-64 49018 143773 30591 305795 189901
# 65+   45243 122401 30830 203634 148843
```

- `age_breaks` are the lower bounds of the groups you want. They must line up with the five-year bands, and the first must be 0. Leaving it out keeps all 20 bands.
- `country` and `regions` select what to include. The order of `regions` is the column order, and the same order is used by `region_distances()`.
- `by_region = FALSE` pools each country into one column, for a model without regions.

The settlement profiles mean age composition differs between regions, so which regions you pick changes the epidemic, not just its size.

## Maps

`read_geography(layer)` reads one of three layers as an `sf` object (needs the `sf` package). The coordinate system is EPSG:32632, with distances in metres.

| Layer | Rows | Contents |
|---|---|---|
| `regions` | 20 | region polygons, with the demographic summary as attributes |
| `settlements` | 540 | points: `name`, `role` (capital, regional centre, town, village), `pop` |
| `cells` | 11,978 | population grid: `pop`, `area_km2`, `density_per_km2`, `settlement` |

`plot_region_map(values)` fills the region polygons by a named vector or a `data.frame(region, value)`, e.g. model output by region. `values = NULL` maps population. Country borders are drawn in black with the country name in the middle, and each region carries its code; turn either off with `country_labels = FALSE` or `region_labels = FALSE`. Regions you give no value for are drawn grey, so a map of one country still shows the others as context.

## Building a region mixing matrix

`build_model()` needs `M_region`, where `M_region[r, s]` is the share of a region-`r` resident's contacts made in region `s`, with rows summing to 1. There's no real travel data here, so you have to assume something. `region_distances()` gives straight-line distances (km) between region centres, which is enough for a simple distance-decay assumption:

```r
dist <- region_distances(country = "Sweway")
within <- 0.9                   # share of contacts made at home
w <- exp(-dist / 100)           # 100 km distance scale
diag(w) <- 0
M <- (1 - within) * w / rowSums(w)
diag(M) <- within
```

The two numbers to think about are `within` (how connected the regions are) and the distance scale. `diag(5)` gives completely isolated regions, which is a useful comparison. See `examples/06_country_data.R`.

## Generating exercise data

`simulate_line_list()` runs an epidemic in one of the countries and returns a line list of cases — the data an exercise hands out:

```r
ll <- simulate_line_list("Sweway", R0 = 1.8, p_detect = 0.1,
                         age_breaks = c(0, 20, 45, 65), seed = 1)
#    id date_onset age_group region vaccinated hospitalised date_hospital  died date_death
# 1:  1 2026-01-01     45-64     S1      FALSE         TRUE    2026-01-07 FALSE       <NA>
# 2:  2 2026-01-01      0-19     S2      FALSE         TRUE    2026-01-23 FALSE       <NA>
# 3:  3 2026-01-04     20-44     S1      FALSE        FALSE          <NA> FALSE       <NA>
```

- The model runs at **full detail** — all 20 five-year bands and all five regions. `age_breaks` only changes how ages are *reported*.
- Each case carries `age_mid`, the middle of its reported age group, so the line list can be re-aggregated into any coarser grouping without parsing labels: `cut(age_mid, c(0, 20, 65, Inf))`. An open-ended top group ("95+", "65+") gets its lower bound plus 2.5.
- **Every hospitalisation and every death is in the list**; other infections appear with probability `p_detect`.
- `R0 = NULL` draws one at random from `R0_range`, so nobody needs to know it in advance. The truth is attached as `attr(ll, "truth")`: the R0 used, `p_detect`, and the model itself.
- The disease is described by arguments, not hard-coded: `latent_period`, `infectious_period`, `susceptibility`, `p_hosp`, `p_death_hosp`, `p_death_comm`, `dur_hosp`, `vacc_coverage` and the three VEs. The defaults are a fast respiratory infection (latent 1 day, infectious 2) with severity rising steeply with age, in an **unvaccinated** population (`vacc_coverage = 0`).
- Each age-varying parameter takes a single number, one value per five-year band, or a **function of age**, so a whole profile fits on one line:

  ```r
  simulate_line_list("Sweway",
    latent_period = 3, infectious_period = 5,
    p_hosp = function(age) pmin(0.002 * exp(0.045 * (age - 10)), 0.5),
    vacc_coverage = function(age) ifelse(age >= 65, 0.8, 0))
  ```

  Anything else (`transmissibility`, `beta_schedule`, …) still goes through `...` to `build_model()`.
- `n_seed`, `seed_region` and `seed_ages` control the seeding; `within` and `scale_km` the mixing between regions. The mixing comes from `region_distances()`, so `simulate_line_list()` needs `sf`.
- If the seeded epidemic fades out, the function retries with a new draw (`max_tries`), and warns if it never takes off.

**Dates** are reconstructed from the model's own timings, each drawn per case: infection to onset is the latent period, onset to admission (or to a community death) is the infectious period, and hospital stays average `dur_hosp`. Hospital deaths are drawn per admission, so they follow the model's death *rates* but not its daily death counts exactly.

`admissions_from_line_list()` turns the line list back into the daily admissions `fit_model()` takes, filling in days with no admissions:

```r
data <- admissions_from_line_list(ll, start_date = min(ll$date_onset))
fit <- fit_model(model, data[time <= 60], pars = "R0", method = "mcmc")
```

`examples/07_line_list_exercise.R` runs this end to end: generate the data, look at it, fit R0 from the first 60 days, and compare with the truth. With default settings the fit recovers R0 to within a few percent.

**A caveat worth repeating to participants:** fitting the early epidemic gives a good estimate of R0, but the *level* of admissions depends on the seeding and the severity assumption together, so `hosp_scale` is poorly identified until the data covers more of the epidemic.
