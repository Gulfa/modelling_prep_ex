# Generating a synthetic line list of cases for exercises.
#
# The model is compartmental, so it gives daily counts, not individuals. This
# runs the stochastic model at full detail (all five-year age bands, all
# regions of a country), then turns the daily event counts into individual
# records: every hospital admission and every death is found, plus a fraction
# of the remaining infections.

#' Simulate a line list of cases
#'
#' Runs an epidemic in one of the bundled countries ([read_demography()]) and
#' returns a line list: one row per detected case, with dates of onset,
#' hospital admission and death. All hospitalised cases and all deaths are
#' detected; other infections are detected with probability `p_detect`.
#'
#' The model runs at full detail (all five-year age bands and all regions of
#' the country); `age_breaks` only controls how ages are reported.
#'
#' The epidemiological parameters are arguments, with defaults describing a
#' respiratory infection whose severity rises steeply with age. Each
#' age-varying one takes a single number, one value per five-year band, or a
#' function of age (evaluated at the midpoint of each band), so a whole
#' profile can be changed in one line:
#'
#' ```
#' simulate_line_list(p_hosp = function(age) 0.01 * exp(0.05 * age),
#'                    vacc_coverage = function(age) pmin(0.9, pmax(0, (age - 20) / 80)))
#' ```
#'
#' Anything else, e.g. `transmissibility` or `beta_schedule`, can still be
#' passed through `...` to [build_model()].
#'
#' Dates are reconstructed from the model's timings: time from infection to
#' onset is the latent period, onset to admission or to a community death is
#' the infectious period, and hospital stays last `dur_hosp` on average, each
#' drawn per case. Hospital deaths are drawn per admission rather than taken
#' from the model's death counts, so deaths match the model's rates but not
#' its daily death counts exactly.
#'
#' The true parameters are attached as `attr(x, "truth")`, so an exercise can
#' hand out the line list and keep the answers.
#'
#' @param country One of the countries in the bundled data.
#' @param R0 Basic reproduction number. `NULL` draws one uniformly from
#'   `R0_range`.
#' @param R0_range Range to draw `R0` from when it is not given.
#' @param t_max Days to simulate.
#' @param p_detect Probability that an infection which is neither hospitalised
#'   nor fatal appears in the line list.
#' @param age_breaks Age groups to report, as in [population_matrix()]
#'   (`NULL` reports the five-year bands the model runs with).
#' @param start_date Date of day 0. `NULL` reports whole days instead of
#'   dates.
#' @param latent_period,infectious_period Mean duration (days) of E and I.
#' @param susceptibility Relative susceptibility by age.
#' @param p_hosp Probability an infection is hospitalised, by age.
#' @param p_death_hosp Probability a hospitalised case dies, by age.
#' @param p_death_comm Probability an infection that is not hospitalised dies,
#'   by age.
#' @param dur_hosp Mean length of hospital stay (days), by age.
#' @param vacc_coverage Fraction vaccinated at the start, by age. The default,
#'   0, gives an unvaccinated population.
#' @param ve_infection,ve_hosp,ve_death Vaccine effectiveness, by age; see
#'   [build_model()].
#' @param n_seed Number of infections seeded at day 0.
#' @param seed_region Region to seed in; `NULL` picks one at random.
#' @param seed_ages Age range the seeded infections are spread over.
#' @param within Share of contacts made in a person's own region; the rest are
#'   spread over the other regions, decaying with distance (`scale_km`).
#' @param scale_km Distance scale (km) for contacts outside the home region.
#' @param seed Random seed.
#' @param dt Time step of the stochastic model.
#' @param min_cases If the seeded epidemic fades out, retry with a new seed
#'   until at least this many cases are detected (at most `max_tries` times).
#' @param max_tries Maximum number of attempts.
#' @param ... Passed to [build_model()], overriding the built-in defaults.
#'
#' @return A `data.table` with one row per detected case: `id`, `date_onset`,
#'   `age_group`, `age_mid`, `region`, `vaccinated`, `hospitalised`,
#'   `date_hospital`, `died` and `date_death`. `age_mid` is the middle of the
#'   reported age group, which makes it easy to re-aggregate into other age
#'   groupings, e.g. `cut(age_mid, c(0, 20, 65, Inf))`. Use
#'   [admissions_from_line_list()] to get the daily admissions that
#'   [fit_model()] takes.
#' @export
#' @examples
#' \dontrun{
#' ll <- simulate_line_list("Sweway", R0 = 1.8, age_breaks = c(0, 20, 45, 65),
#'                          seed = 1)
#' head(ll)
#' attr(ll, "truth")$R0
#'
#' # a slower disease, with the elderly vaccinated
#' ll2 <- simulate_line_list(
#'   "Denden", R0 = 2, latent_period = 3, infectious_period = 5,
#'   vacc_coverage = function(age) ifelse(age >= 65, 0.8, 0),
#'   seed = 1)
#' }
simulate_line_list <- function(country = "Sweway", R0 = NULL,
                               R0_range = c(1.3, 2.5), t_max = 180,
                               p_detect = 0.1, age_breaks = NULL,
                               start_date = as.Date("2026-01-01"),
                               latent_period = 1, infectious_period = 2,
                               susceptibility = function(age)
                                 ifelse(age < 20, 0.7, 1),
                               p_hosp = function(age)
                                 pmin(0.004 * exp(0.045 * (age - 10)), 0.5),
                               p_death_hosp = function(age)
                                 pmin(0.1 * exp(0.030 * (age - 10)), 0.5),
                               p_death_comm = function(age)
                                 pmin(0.0004 * exp(0.030 * (age - 10)), 0.05),
                               dur_hosp = 7,
                               vacc_coverage = 0,
                               ve_infection = 0.4, ve_hosp = 0.7,
                               ve_death = 0.5,
                               n_seed = 20, seed_region = NULL,
                               seed_ages = c(20, 60),
                               within = 0.9, scale_km = 100, seed = NULL,
                               dt = 0.25, min_cases = 50, max_tries = 10,
                               ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  check_range(p_detect, "p_detect", max = 1)
  if (is.null(R0)) {
    R0 <- stats::runif(1, R0_range[1], R0_range[2])
  }
  epi <- list(latent_period = latent_period,
              infectious_period = infectious_period,
              susceptibility = susceptibility, p_hosp = p_hosp,
              p_death_hosp = p_death_hosp, p_death_comm = p_death_comm,
              dur_hosp = dur_hosp, vacc_coverage = vacc_coverage,
              ve_infection = ve_infection, ve_hosp = ve_hosp,
              ve_death = ve_death)
  model <- line_list_model(country, R0, epi, n_seed, seed_region, seed_ages,
                           within, scale_km, ...)

  for (try in seq_len(max_tries)) {
    ll <- line_list_from_model(model, t_max, p_detect, dt, start_date)
    if (nrow(ll) >= min_cases) {
      break
    }
  }
  if (nrow(ll) < min_cases) {
    warning(sprintf(paste("The epidemic faded out in all %d attempts; only %d",
                          "case(s) detected. Try a larger R0 or n_seed."),
                    max_tries, nrow(ll)))
  }
  if (!is.null(age_breaks)) {
    ll[, age_group := regroup_ages(age_group, age_breaks)]
  }
  # middle of the reported age group, so cases can be re-aggregated into any
  # coarser grouping, e.g. cut(age_mid, c(0, 20, 65, Inf))
  ll[, age_mid := age_middle(age_group)]
  setcolorder(ll, c("id", "date_onset", "age_group", "age_mid"))
  setattr(ll, "truth",
          list(country = country, R0 = R0, p_detect = p_detect,
               start_date = start_date, model = model))
  ll[]
}

# The model used by simulate_line_list(): full age and region detail for one
# country. `epi` holds the arguments of simulate_line_list() that describe the
# disease; the age-varying ones may be given as a function of age.
line_list_model <- function(country, R0, epi, n_seed, seed_region, seed_ages,
                            within, scale_km, ...) {
  N <- population_matrix(country)
  ages <- rownames(N)
  mid <- age_middle(ages)

  if (is.null(seed_region)) {
    seed_region <- sample(colnames(N), 1)
  } else if (!seed_region %in% colnames(N)) {
    stop(sprintf("'seed_region' must be a region of %s: %s", country,
                 paste(colnames(N), collapse = ", ")))
  }
  if (length(seed_ages) != 2 || seed_ages[1] >= seed_ages[2]) {
    stop("'seed_ages' must be two ages, youngest first")
  }
  seeded <- which(mid > seed_ages[1] & mid < seed_ages[2])
  if (length(seeded) == 0) {
    stop("No age band falls inside 'seed_ages'")
  }
  I0 <- matrix(0, nrow(N), ncol(N), dimnames = dimnames(N))
  I0[seeded, seed_region] <- stats::rmultinom(1, n_seed, rep(1, length(seeded)))

  by_age <- c("susceptibility", "p_hosp", "p_death_hosp", "p_death_comm",
              "dur_hosp", "vacc_coverage", "ve_infection", "ve_hosp",
              "ve_death")
  epi[by_age] <- lapply(by_age, function(nm) age_param(epi[[nm]], mid, nm))

  args <- c(list(N = N,
                 C_age = age_contact_matrix(ages),
                 M_region = distance_mixing(country, within, scale_km),
                 R0 = R0,
                 I0 = I0),
            epi)
  args[names(list(...))] <- list(...)
  do.call(build_model, args)
}

# An age-varying parameter: a function of age, one value per age band, or a
# single value for all ages.
age_param <- function(x, mid, name) {
  if (is.function(x)) {
    x <- x(mid)
  }
  if (!is.numeric(x) || anyNA(x) || !length(x) %in% c(1, length(mid))) {
    stop(sprintf(paste("'%s' must be a number, a vector with one value per",
                       "five-year age band (%d), or a function of age"),
                 name, length(mid)))
  }
  x
}

#' A synthetic age contact matrix
#'
#' Contacts per day between age groups, made up of contacts with people of
#' similar age plus a parent/child component one generation apart. It is a
#' stand-in for a real contact survey; only its shape matters when
#' [build_model()] is given `R0`, not its scale.
#'
#' @param ages Age group labels, e.g. `rownames(population_matrix(...))`.
#' @param width Width (years) of the similar-age component.
#' @param generation Age difference (years) of the parent/child component.
#' @param generation_weight Size of the parent/child component relative to the
#'   similar-age one.
#'
#' @return A square matrix with `ages` as dimnames.
#' @export
#' @examples
#' round(age_contact_matrix(c("0-19", "20-44", "45-64", "65+")), 2)
age_contact_matrix <- function(ages, width = 8, generation = 30,
                               generation_weight = 0.4) {
  mid <- age_lower(ages)
  mid <- mid + c(diff(mid), 5) / 2
  d <- outer(mid, mid, "-")
  m <- exp(-(d / width)^2) +
    generation_weight * (exp(-((d - generation) / width)^2) +
                           exp(-((d + generation) / width)^2))
  dimnames(m) <- list(ages, ages)
  m
}

# Region mixing from the distance between region centres: `within` of contacts
# at home, the rest decaying with distance.
distance_mixing <- function(country, within = 0.9, scale_km = 100) {
  d <- region_distances(country = country)
  if (nrow(d) == 1) {
    return(matrix(1, 1, 1, dimnames = dimnames(d)))
  }
  check_range(within, "within", max = 1)
  w <- exp(-d / scale_km)
  diag(w) <- 0
  m <- (1 - within) * w / rowSums(w)
  diag(m) <- within
  m
}

# Turn one stochastic run into individual case records
line_list_from_model <- function(model, t_max, p_detect, dt, start_date) {
  res <- run_model(model, "stochastic", t_max = t_max, n_sim = 1, dt = dt,
                   vars = c("inc_infections", "inc_hosp", "inc_deaths_comm"))
  d <- get_output(res, by = c("age", "region", "vacc_status"))
  d <- dcast(d, time + age + region + vacc_status ~ variable,
             value.var = "value")
  d <- d[inc_infections + inc_hosp + inc_deaths_comm > 0]

  # per-group probabilities, in the backend's group order
  pars <- model$pars
  group <- rep(seq_along(model$N_group), 2)
  vacc <- rep(1:2, each = length(model$N_group))
  key <- paste(rep(rep(model$dims$age, length(model$dims$region)), 2),
               rep(rep(model$dims$region, each = length(model$dims$age)), 2),
               c("unvaccinated", "vaccinated")[vacc])
  idx <- match(paste(d$age, d$region, d$vacc_status), key)
  p_hosp <- pars$p_hosp[cbind(group, vacc)][idx]
  p_death_comm <- pars$p_death_comm[cbind(group, vacc)][idx]
  p_death_hosp <- pars$p_death_hosp[cbind(group, vacc)][idx]
  dur_hosp <- pars$dur_hosp[cbind(group, vacc)][idx]

  # Mild cases are detected among infections that do not become severe, so
  # they are not counted twice with the hospitalisations and deaths below.
  p_severe <- p_hosp + (1 - p_hosp) * p_death_comm
  n_mild <- stats::rbinom(nrow(d), as.integer(d$inc_infections),
                          p_detect * (1 - p_severe))

  latent <- 1 / pars$sigma
  infectious <- 1 / pars$gamma
  rows <- rbind(
    # detected infections that stay mild: onset is after the latent period
    expand_cases(d, n_mild, kind = "mild", latent = latent),
    expand_cases(d, as.integer(d$inc_hosp), kind = "hospital",
                 infectious = infectious, p_death = p_death_hosp,
                 dur_hosp = dur_hosp),
    expand_cases(d, as.integer(d$inc_deaths_comm), kind = "community_death",
                 infectious = infectious))
  if (is.null(rows) || nrow(rows) == 0) {   # the epidemic faded out
    return(empty_line_list(start_date))
  }

  setorder(rows, date_onset)
  rows[, id := seq_len(.N)]
  if (!is.null(start_date)) {
    for (v in c("date_onset", "date_hospital", "date_death")) {
      set(rows, j = v, value = start_date + rows[[v]])
    }
  }
  setcolorder(rows, c("id", "date_onset", "age_group", "region", "vaccinated",
                      "hospitalised", "date_hospital", "died", "date_death"))
  rows[]
}

# One row per case, with dates drawn from the model's delay distributions
expand_cases <- function(d, n, kind, latent = NULL, infectious = NULL,
                         p_death = NULL, dur_hosp = NULL) {
  keep <- which(n > 0)
  if (length(keep) == 0) {
    return(NULL)
  }
  i <- rep(keep, n[keep])
  total <- length(i)
  out <- data.table(
    age_group = d$age[i],
    region = d$region[i],
    vaccinated = d$vacc_status[i] == "vaccinated",
    hospitalised = kind == "hospital",
    died = FALSE,
    date_hospital = NA_real_,
    date_death = NA_real_)

  if (kind == "mild") {
    # infected on this day, onset after the latent period
    out[, date_onset := d$time[i] + round(stats::rexp(total, 1 / latent))]
  } else {
    # hospitalised, or died in the community, at the end of the infectious
    # period: work backwards to onset
    out[, date_onset := pmax(0, d$time[i] -
                               round(stats::rexp(total, 1 / infectious)))]
    if (kind == "hospital") {
      out[, date_hospital := d$time[i]]
      has_died <- stats::rbinom(total, 1, p_death[i]) == 1
      stay <- round(stats::rexp(total, 1 / dur_hosp[i]))
      out[, died := has_died]
      out[has_died, date_death := date_hospital + stay[has_died]]
    } else {
      out[, `:=`(died = TRUE, date_death = d$time[i])]
    }
  }
  out
}

empty_line_list <- function(start_date) {
  d <- if (is.null(start_date)) numeric() else start_date[0]
  data.table(id = integer(), date_onset = d, age_group = character(),
             age_mid = numeric(), region = character(), vaccinated = logical(),
             hospitalised = logical(), date_hospital = d, died = logical(),
             date_death = d)
}

#' Daily hospital admissions from a line list
#'
#' Counts admissions per day, in the format [fit_model()] expects: a `time`
#' column counting days from the start of the data, and `admissions`.
#'
#' @param line_list A line list from [simulate_line_list()].
#' @param by Also split the series by `"region"` and/or `"age_group"`.
#' @param start_date Day 1 of the output. The default is the first admission.
#'   Pass the date the model should start from if you want the epidemic to be
#'   already running when the data begins.
#'
#' @return A `data.table` with `time`, the `by` columns and `admissions`,
#'   including days with no admissions.
#' @export
#' @examples
#' \dontrun{
#' ll <- simulate_line_list("Sweway", R0 = 1.8, seed = 1)
#' data <- admissions_from_line_list(ll)
#' fit_model(model, data, pars = "R0")
#' }
admissions_from_line_list <- function(line_list, by = character(0),
                                      start_date = NULL) {
  bad <- setdiff(by, c("region", "age_group"))
  if (length(bad) > 0) {
    stop("'by' must be a subset of: region, age_group")
  }
  d <- as.data.table(line_list)[!is.na(date_hospital)]
  if (nrow(d) == 0) {
    stop("The line list has no hospital admissions")
  }
  if (is.null(start_date)) {
    start_date <- min(d$date_hospital)
  } else if (length(start_date) != 1 || is.na(start_date)) {
    stop("'start_date' must be a single, non-missing date")
  }
  d[, time := as.numeric(date_hospital - start_date) + 1]
  if (any(d$time < 1)) {
    stop("'start_date' is after the first admission in the line list")
  }
  counts <- d[, .(admissions = .N), by = c("time", by)]
  # fill in the days with no admissions
  grid <- do.call(CJ, c(list(time = seq_len(max(counts$time))),
                        lapply(by, function(v) sort(unique(d[[v]])))))
  if (length(by) > 0) {
    setnames(grid, c("time", by))
  }
  out <- counts[grid, on = c("time", by)]
  out[is.na(admissions), admissions := 0L]
  setorderv(out, c(by, "time"))
  out[]
}
