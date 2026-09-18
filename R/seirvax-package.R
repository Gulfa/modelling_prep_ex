#' seirvax: SEIR model with hospitalisation, death and vaccination
#'
#' Build a model with [build_model()], run it with [run_model()], and look at
#' the results with [get_output()], `summary()` and `plot()`. See
#' `docs/model_description.md` for the model equations.
#'
#' @useDynLib seirvax, .registration = TRUE
#' @import data.table
#' @keywords internal
"_PACKAGE"

# data.table non-standard evaluation
utils::globalVariables(c(".", "value", "variable", "time", "sim", "measure",
                         "median", "lower", "upper", ".SD", ".trace"))
