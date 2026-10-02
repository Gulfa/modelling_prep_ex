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
  t_end <- max(object$time)
  final <- get_output(object, c("cum_infections", "cum_hosp", "D"), by)
  final <- final[time == t_end]
  final[, measure := c(cum_infections = "infections", cum_hosp = "hosp_admissions",
                       D = "deaths")[variable]]

  H <- get_output(object, "H", by)
  peak <- H[, .(peak_hosp_occupancy = max(value),
                peak_hosp_day = time[which.max(value)]), by = c("sim", by)]
  peak <- melt(peak, id.vars = c("sim", by), variable.name = "measure",
               variable.factor = FALSE)

  per_sim <- rbind(final[, c("sim", by, "measure", "value"), with = FALSE],
                   peak[, c("sim", by, "measure", "value"), with = FALSE])
  measures <- c("infections", "hosp_admissions", "deaths",
                "peak_hosp_occupancy", "peak_hosp_day")
  per_sim[, measure := factor(measure, measures)]

  if (object$type == "deterministic") {
    ret <- per_sim[, c(by, "measure", "value"), with = FALSE]
  } else {
    a <- (1 - level) / 2
    ret <- per_sim[, .(mean = mean(value), median = stats::median(value),
                       lower = stats::quantile(value, a, names = FALSE),
                       upper = stats::quantile(value, 1 - a, names = FALSE)),
                   by = c(by, "measure")]
  }
  setorderv(ret, c("measure", by))
  ret[, measure := as.character(measure)]
  ret[]
}

#' Plot a model run
#'
#' Plots one output variable over time. Deterministic runs are drawn as
#' lines; stochastic runs as the median with a central interval, optionally
#' with individual simulations.
#'
#' @param x Result of [run_model()].
#' @param var Output variable, see [get_output()].
#' @param by Stratifications, as in [get_output()]. The first is shown by
#'   colour, any others as facets.
#' @param level Width of the interval for stochastic runs.
#' @param n_traces Number of individual stochastic simulations to draw.
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
plot.seirvax_run <- function(x, var = "H", by = character(0), level = 0.95,
                             n_traces = 0, ...) {
  if (length(var) != 1) {
    stop("'var' must be a single variable")
  }
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
    a <- (1 - level) / 2
    q <- d[, .(median = stats::median(value),
               lower = stats::quantile(value, a, names = FALSE),
               upper = stats::quantile(value, 1 - a, names = FALSE)),
           by = c("time", by)]
    p <- ggplot2::ggplot(q, ggplot2::aes(time)) + map_colour
    if (n_traces > 0) {
      sims <- seq_len(min(n_traces, x$n_sim))
      tr <- d[sim %in% sims]
      tr[, .trace := do.call(paste, c(list(sim), .SD)), .SDcols = by]
      p <- p + ggplot2::geom_line(ggplot2::aes(y = value, group = .trace),
                                  data = tr, alpha = 0.2, linewidth = 0.3)
    }
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = 0.25, colour = NA) +
      ggplot2::geom_line(ggplot2::aes(y = median))
  }
  if (length(by) > 1) {
    p <- p + ggplot2::facet_wrap(by[-1], scales = "free_y")
  }
  p + ggplot2::labs(x = "Day", y = var) + ggplot2::theme_bw()
}
