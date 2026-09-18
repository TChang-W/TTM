get_weights <- function(dat, dropout_formula){

  data_death <- dat[dat$death == 1, , drop = FALSE]
  censor_OS <- sort(unique(dat$OS_time[dat$event == 0]))

  dropout_model <- survival::coxph(dropout_formula, data = dat, x = TRUE, y = TRUE)
  drop_basehaz <- survival::basehaz(dropout_model, centered = FALSE)

  id_vec <- unique(data_death$id)
  drop_cumu_haz_vec <- numeric(length(id_vec))
  drop_weight_vec <- numeric(length(id_vec))
  censor_weight_vec <- numeric(length(id_vec))

  for (j in seq_along(id_vec)) {
    i <- id_vec[j]
    data_pt <- data_death[data_death$id == i, , drop = FALSE]
    data_pt <- data_pt[1, , drop = FALSE]

    death_time_i <- data_pt$OS_time
    lp_i <- stats::predict(dropout_model, newdata = data_pt, type = "lp")
    risk_i <- exp(lp_i)

    H0_idx <- drop_basehaz$time <= death_time_i
    H0_ti <- if (any(H0_idx)) max(drop_basehaz$hazard[H0_idx], na.rm = TRUE) else 0
    if (!is.finite(H0_ti)) H0_ti <- 0

    drop_cumu_haz_vec[j] <- H0_ti * risk_i
    drop_weight_vec[j] <- exp(-drop_cumu_haz_vec[j])
    censor_weight_vec[j] <- mean(censor_OS >= data_pt$OS_time)
  }

  weights_df <- data.frame(
    id = id_vec,
    weights = drop_weight_vec * censor_weight_vec,
    stringsAsFactors = FALSE
  )
  return(weights_df)
}
