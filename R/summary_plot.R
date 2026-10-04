#' Summarise a model run
#'
#' Key outcomes per simulation: total infections, hospital admissions and
#' deaths by the end of the run, and the peak hospital occupancy and its
#' timing. For stochastic runs these are summarised over simulations by the
#' mean, median and a central interval.
#'
#' @param object Result of [run_model()].
#' @param by Stratifications, as in [get_output()].
#' @param level Width of the interval for stochastic runs.
#' @param ... Unused.
#'
#' @return A `data.table` with one row per outcome (and `by` group).
#' @export
summary.seirvax_run <- function(object, by = character(0), level = 0.95, ...) {
  by <- check_by(by)
  d <- get_output(object, c("cum_infections", "cum_hosp", "D", "H"), by)
  per_sim <- trajectory_measures(d, "sim", by)

  if (object$type == "deterministic") {
    ret <- per_sim[, c(by, "measure", "value"), with = FALSE]
  } else {
    ret <- summarise_measures(per_sim, by, level)
  }
  setorderv(ret, c("measure", by))
  ret[, measure := as.character(measure)]
  ret[]
}

outcome_measures <- c("infections", "hosp_admissions", "deaths",
                      "peak_hosp_occupancy", "peak_hosp_day")

# Outcomes of each trajectory (identified by the `id` columns): totals at the
# end of the run, and the peak of hospital occupancy with its timing.
trajectory_measures <- function(d, id, by) {
  final <- d[variable %in% c("cum_infections", "cum_hosp", "D") &
               time == max(time)]
  final[, measure := c(cum_infections = "infections",
                       cum_hosp = "hosp_admissions", D = "deaths")[variable]]
  peak <- d[variable == "H",
            .(peak_hosp_occupancy = max(value),
              peak_hosp_day = time[which.max(value)]), by = c(id, by)]
  peak <- melt(peak, id.vars = c(id, by), variable.name = "measure",
               variable.factor = FALSE)
  cols <- c(id, by, "measure", "value")
  ret <- rbind(final[, cols, with = FALSE], peak[, cols, with = FALSE])
  ret[, measure := factor(measure, outcome_measures)]
  ret[]
}

# Median, mean and a central interval over trajectories
summarise_measures <- function(per_traj, by, level) {
  a <- (1 - level) / 2
  per_traj[, .(mean = mean(value), median = stats::median(value),
               lower = stats::quantile(value, a, names = FALSE),
               upper = stats::quantile(value, 1 - a, names = FALSE)),
           by = c(by, "measure")]
}

#' Plot a model run
#'
#' Plots one output variable over time. Deterministic runs are drawn as
#' lines. Stochastic runs are drawn as a spaghetti plot by default: one line
#' per simulation, which shows how variable the epidemics are, including ones
#' that fade out or peak at a very different time. `style = "ribbon"` shows
#' the median and a central interval instead, which is tidier but hides that
#' individual epidemics never look like the median.
#'
#' @param x Result of [run_model()].
#' @param var Output variable, see [get_output()].
#' @param by Stratifications, as in [get_output()]. The first is shown by
#'   colour, any others as facets.
#' @param style `"spaghetti"` (one line per simulation) or `"ribbon"`
#'   (median and interval). Ignored for deterministic runs.
#' @param n_traces Maximum number of simulations to draw as individual lines.
#'   The default is 50 for `style = "spaghetti"` and none for `"ribbon"`,
#'   where you can still ask for some.
#' @param show_median Draw the median on top of the spaghetti.
#' @param level Width of the interval for `style = "ribbon"`.
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' m <- build_model(N = 1e4, R0 = 1.6, I0 = 5, p_hosp = 0.05)
#' res <- run_model(m, "stochastic", t_max = 150, n_sim = 100, seed = 1)
#' plot(res, var = "I")
#' plot(res, var = "I", style = "ribbon")
plot.seirvax_run <- function(x, var = "H", by = character(0),
                             style = c("spaghetti", "ribbon"),
                             n_traces = NULL, show_median = TRUE,
                             level = 0.95, ...) {
  if (length(var) != 1) {
    stop("'var' must be a single variable")
  }
  style <- match.arg(style)
  n_traces <- n_traces %||% if (style == "spaghetti") 50 else 0
  by <- check_by(by)
  d <- get_output(x, var, by)
  map_colour <- if (length(by) > 0) {
    ggplot2::aes(colour = .data[[by[[1]]]], fill = .data[[by[[1]]]])
  } else {
    ggplot2::aes()
  }

  if (x$type == "deterministic") {
    p <- ggplot2::ggplot(d, ggplot2::aes(time, value)) + map_colour +
      ggplot2::geom_line()
  } else {
    traces <- function(p, alpha) {
      tr <- d[sim %in% seq_len(min(n_traces, x$n_sim))]
      tr[, .trace := do.call(paste, c(list(sim), .SD)), .SDcols = by]
      p + ggplot2::geom_line(ggplot2::aes(y = value, group = .trace),
                             data = tr, alpha = alpha, linewidth = 0.3)
    }
    q <- NULL
    if (style == "ribbon" || show_median) {
      a <- (1 - level) / 2
      q <- d[, .(median = stats::median(value),
                 lower = stats::quantile(value, a, names = FALSE),
                 upper = stats::quantile(value, 1 - a, names = FALSE)),
             by = c("time", by)]
    }
    p <- ggplot2::ggplot(q %||% d, ggplot2::aes(time)) + map_colour
    if (style == "spaghetti") {
      # fainter lines the more there are, so the bulk stays readable
      p <- traces(p, max(0.08, min(0.6, 15 / min(n_traces, x$n_sim))))
      if (show_median) {
        p <- p + ggplot2::geom_line(ggplot2::aes(y = median), linewidth = 0.8)
      }
    } else {
      if (n_traces > 0) {
        p <- traces(p, 0.2)
      }
      p <- p +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                             alpha = 0.25, colour = NA) +
        ggplot2::geom_line(ggplot2::aes(y = median))
    }
  }
  if (length(by) > 1) {
    p <- p + ggplot2::facet_wrap(by[-1], scales = "free_y")
  }
  p + ggplot2::labs(x = "Day", y = var) + ggplot2::theme_bw()
}
