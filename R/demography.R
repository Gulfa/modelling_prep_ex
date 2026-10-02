# Reading and selecting the bundled demographic and geographic data for the
# four exercise countries (Sweway, Finmark, Norland, Denden).

#' Demographic data for the exercise countries
#'
#' Population by country, region and five-year age group, with region area,
#' density and settlement profile. The data are fictional, made for
#' preparedness exercises.
#'
#' @param file Path to the data; defaults to the copy bundled with the
#'   package.
#'
#' @return A `data.table` with one row per region and age group:
#'   `country`, `region`, `age_group` (an ordered factor), `population`,
#'   `share_of_region`, `area_km2`, `density_per_km2`, `pop_total` and
#'   `profile`.
#' @export
#' @examples
#' d <- read_demography()
#' head(d)
#' region_info()
read_demography <- function(file = NULL) {
  if (is.null(file)) {
    file <- system.file("extdata", "region_demographics_model_long.csv",
                        package = "seirvax")
  }
  if (!nzchar(file) || !file.exists(file)) {
    stop("Demographic data not found; give 'file'")
  }
  d <- fread(file)
  setnames(d, tolower(names(d)))
  d[, age_group := factor(age_group, unique(age_group[order(age_lower(age_group))]),
                          ordered = TRUE)]
  setkey(d, country, region)
  d[]
}

# Lower bound of an age group label: "0-4" -> 0, "95+" -> 95
age_lower <- function(x) {
  as.numeric(sub("[-+].*$", "", as.character(x)))
}

#' One row per region
#'
#' @param country Countries to keep (default: all).
#' @param regions Regions to keep (default: all regions of `country`). The
#'   order given here is the order used by [population_matrix()] and
#'   [region_distances()].
#' @param demography Data from [read_demography()]; read if not given.
#'
#' @return A `data.table` with `country`, `region`, `profile`, `area_km2`,
#'   `density_per_km2` and `pop_total`.
#' @export
region_info <- function(country = NULL, regions = NULL, demography = NULL) {
  d <- demography %||% read_demography()
  d <- select_regions(d, country, regions)
  unique(d[, .(country, region, profile, area_km2, density_per_km2, pop_total)])
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Filter to the requested countries / regions, keeping the requested order
# (and otherwise the order in the data file).
select_regions <- function(d, country = NULL, regions = NULL) {
  d <- copy(d)   # never modify the caller's data by reference
  keep_country <- country   # not 'country': that is also a column name
  if (!is.null(country)) {
    unknown <- setdiff(country, unique(d$country))
    if (length(unknown) > 0) {
      stop(sprintf("Unknown country: %s. Available: %s",
                   paste(unknown, collapse = ", "),
                   paste(unique(d$country), collapse = ", ")))
    }
    d <- d[country %chin% keep_country]
  }
  if (!is.null(regions)) {
    unknown <- setdiff(regions, unique(d$region))
    if (length(unknown) > 0) {
      stop(sprintf("Unknown region: %s. Available: %s",
                   paste(unknown, collapse = ", "),
                   paste(unique(d$region), collapse = ", ")))
    }
    d <- d[region %chin% regions]
    d[, region := factor(region, regions)]
    setorder(d, region)
    d[, region := as.character(region)]
  }
  if (nrow(d) == 0) {
    stop("No regions selected")
  }
  d
}

#' Population matrix for a model
#'
#' Builds the age × region population matrix that [build_model()] takes as
#' `N`, for a choice of countries, regions and age grouping.
#'
#' @inheritParams region_info
#' @param age_breaks Lower bounds of the age groups to aggregate the
#'   five-year bands into, e.g. `c(0, 20, 45, 65)` gives `0-19`, `20-44`,
#'   `45-64` and `65+`. Must be multiples of 5 and at most 95. `NULL` (the
#'   default) keeps all five-year bands; a single value, such as `0`, gives
#'   one age group.
#' @param by_region If `FALSE`, regions are pooled per country, giving one
#'   column per country.
#'
#' @return A matrix of populations, age groups in rows and regions (or
#'   countries) in columns, with dimnames.
#' @export
#' @examples
#' # One country, four age groups, five regions
#' N <- population_matrix("Sweway", age_breaks = c(0, 20, 45, 65))
#' N
#'
#' # Two regions, all five-year bands
#' population_matrix(regions = c("S1", "S2"))
#'
#' # Whole countries as single populations
#' population_matrix(by_region = FALSE, age_breaks = c(0, 65))
population_matrix <- function(country = NULL, regions = NULL,
                              age_breaks = NULL, by_region = TRUE,
                              demography = NULL) {
  d <- demography %||% read_demography()
  d <- select_regions(d, country, regions)
  d[, age_group := regroup_ages(age_group, age_breaks)]
  area <- if (by_region) "region" else "country"
  agg <- d[, .(population = sum(population)),
           by = c("age_group", area)]
  cols <- unique(d[[area]])
  m <- matrix(0, nlevels(agg$age_group), length(cols),
              dimnames = list(levels(agg$age_group), cols))
  m[cbind(as.integer(agg$age_group), match(agg[[area]], cols))] <- agg$population
  m
}

# Collapse five-year age bands into groups starting at `breaks`
regroup_ages <- function(age_group, breaks) {
  levs <- levels(age_group)
  if (is.null(breaks)) {
    return(age_group)
  }
  breaks <- sort(unique(breaks))
  lower <- age_lower(levs)
  if (any(!breaks %in% lower)) {
    stop(sprintf(paste("'age_breaks' must be the start of an age band in the",
                       "data: multiples of 5 up to %d"), max(lower)))
  }
  idx <- findInterval(lower, breaks)
  if (any(idx == 0)) {
    stop("The first value of 'age_breaks' must be 0")
  }
  labels <- vapply(seq_along(breaks), function(k) {
    if (k == length(breaks)) {
      sprintf("%d+", breaks[k])
    } else {
      sprintf("%d-%d", breaks[k], breaks[k + 1] - 1)
    }
  }, "")
  # idx is per level of age_group; map it through to each element
  factor(labels[idx][as.integer(age_group)], labels, ordered = TRUE)
}

#' Map data for the exercise countries
#'
#' Reads a layer of the bundled GeoPackage. Needs the `sf` package.
#'
#' @inheritParams region_info
#' @param layer `"regions"` (polygons, with the region's demographics),
#'   `"settlements"` (points, with name, role and population) or `"cells"`
#'   (points on a population grid).
#' @param file Path to the GeoPackage; defaults to the bundled copy.
#'
#' @return An `sf` data.frame, filtered to the selected countries/regions.
#' @export
#' @examples
#' \dontrun{
#' regions <- read_geography("regions", country = "Sweway")
#' plot(regions["POP_TOTAL"])
#' }
read_geography <- function(layer = c("regions", "settlements", "cells"),
                           country = NULL, regions = NULL, file = NULL) {
  layer <- match.arg(layer)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is needed to read the map data")
  }
  if (is.null(file)) {
    file <- system.file("extdata", "geography_model.gpkg", package = "seirvax")
  }
  if (!nzchar(file) || !file.exists(file)) {
    stop("Map data not found; give 'file'")
  }
  d <- sf::st_read(file, layer = layer, quiet = TRUE)
  # the layers use upper case names in 'regions', lower case in the others
  names(d)[names(d) %in% c("COUNTRY", "REGION")] <-
    tolower(names(d)[names(d) %in% c("COUNTRY", "REGION")])
  if (!is.null(country)) {
    d <- d[d$country %in% country, ]
  }
  if (!is.null(regions)) {
    d <- d[d$region %in% regions, ]
    d <- d[order(match(d$region, regions)), ]
  }
  if (nrow(d) == 0) {
    stop("No features selected")
  }
  d
}

#' Distances between region centres
#'
#' Straight-line distances between the centroids of the region polygons, in
#' kilometres. Useful as a starting point for a region mixing matrix (for
#' example a gravity model), which [build_model()] takes as `M_region`.
#'
#' @inheritParams region_info
#' @param file Path to the GeoPackage; defaults to the bundled copy.
#'
#' @return A symmetric matrix with zero diagonal, in the same region order as
#'   [population_matrix()].
#' @export
#' @examples
#' \dontrun{
#' round(region_distances(country = "Sweway"))
#' }
region_distances <- function(country = NULL, regions = NULL, file = NULL) {
  d <- read_geography("regions", country, regions, file)
  centroids <- sf::st_centroid(sf::st_geometry(d))
  m <- sf::st_distance(centroids)
  m <- matrix(as.numeric(m) / 1000, nrow(d), nrow(d),
              dimnames = list(d$region, d$region))
  m
}

#' Map a value by region
#'
#' Fills the region polygons by a value, e.g. model output by region.
#'
#' @param values A named vector (names are regions), or a data.frame with a
#'   `region` column and one other column holding the values. `NULL` fills
#'   by population.
#' @inheritParams region_info
#' @param label Legend title.
#' @param labels Draw the region code on each region.
#' @param file Path to the GeoPackage; defaults to the bundled copy.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' res <- run_model(build_model(N = population_matrix("Sweway"), ...), t_max = 100)
#' d <- get_output(res, "cum_infections", by = "region")[time == 100]
#' plot_region_map(d[, .(region, value)], label = "Infections")
#' }
plot_region_map <- function(values = NULL, country = NULL, regions = NULL,
                            label = NULL, labels = TRUE, file = NULL) {
  d <- read_geography("regions", country, regions, file)
  if (is.null(values)) {
    d$.value <- d$POP_TOTAL
    label <- label %||% "Population"
  } else {
    if (is.data.frame(values)) {
      if (!"region" %in% names(values) || ncol(values) != 2) {
        stop("'values' must be a data.frame with a 'region' column and one ",
             "value column, or a named vector")
      }
      value_col <- setdiff(names(values), "region")
      label <- label %||% value_col
      values <- stats::setNames(values[[value_col]], values$region)
    }
    unknown <- setdiff(names(values), d$region)
    if (is.null(names(values)) || length(unknown) > 0) {
      stop("'values' must be named with regions shown on the map")
    }
    d$.value <- unname(values[d$region])
  }
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data$.value), colour = "white") +
    ggplot2::scale_fill_viridis_c(name = label %||% "Value") +
    ggplot2::theme_void()
  if (labels) {
    p <- p + ggplot2::geom_sf_label(ggplot2::aes(label = .data$region),
                                    size = 3, label.size = 0, alpha = 0.7)
  }
  p
}
