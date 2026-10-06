#' @keywords internal
#' @importFrom stats gaussian median na.omit pnorm quantile
"_PACKAGE"

# Column names used with non-standard evaluation (dplyr / ggplot2).
utils::globalVariables(c(
  "id", "OS_time", "t_pros", "t_retro", "death", "hazard", "death_imputed",
  "time_grid", "estimate", "lower", "upper"
))
