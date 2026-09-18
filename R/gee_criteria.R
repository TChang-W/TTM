gee_criteria <- function(knots_results){
  metrics <- sapply(knots_results, function(kr) {
    geer::geecriteria(kr$gee_model)
  })
  # Transpose metrics
  metrics <- t(metrics)

  metrics_df <- data.frame(
    n_inner_knot = sapply(knots_results, function(kr) kr$n_inner_knot),
    metrics
  )

  return(metrics_df)
}

