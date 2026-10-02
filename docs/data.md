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
