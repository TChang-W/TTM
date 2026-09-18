##' IMPORTANT: This function is for TTM_linear object with B-spline only.
##' Compute the change in slope of the either baseline trajectory or treatment-effect
##' around a knot and test its significance.
##'
##' @param new_times Numeric vector of knot positions. Needs to be of length 3 where the
##'   middle value is the knot and the other two values are slightly before and after the knot.
##' @param fit A fitted `TTM_linear()` object
##' @param var_name   Either `"X0"` (baseline trajectory beta_mu) or
##' `"A"` (treatment effect trajectory beta_A)
##' @return A data.frame with slope1, slope2, diff, se_diff, p_value
##' @export
test_beta_slope <- function(new_times, fit, var_name) {
  res_df <- data.frame()


  results  <- fit

  if (var_name == "A") {
    beta_hat   <- results$weighted_gee$fixed_effects
    beta_vcov  <- results$weighted_gee$fixed_effects_vcov
    var_names <- c("A", grep("^AX_spline", names(beta_hat), value = TRUE))
  } else if (var_name == "X0") {
    beta_hat   <- results$weighted_gee$fixed_effects
    beta_vcov  <- results$weighted_gee$fixed_effects_vcov
    var_names <- c("X0", grep("^X_spline", names(beta_hat), value = TRUE))
  } else {
    stop("var_name must be either \"A\" or \"X0\"")
  }
  beta_coef  <- beta_hat[var_names]

  spline_val <- splines::bs(
    new_times,
    knots   = results$spline_inner_knots,
    degree  = 1,
    intercept = FALSE
  )

  X_mat <- if (var_name == "A") {
    cbind(A = 1, spline_val)
  } else {
    cbind(X0 = 1, spline_val)
  }

  beta_vals <- as.vector(X_mat %*% beta_coef)

  slope1 <- round((beta_vals[2] - beta_vals[1])/(new_times[2] - new_times[1]), 3)
  slope2 <- round((beta_vals[3] - beta_vals[2])/(new_times[3] - new_times[2]), 3)
  diff   <- round(slope2 - slope1, 2)

  x_vec  <- as.matrix(X_mat[1, ] + X_mat[3, ] - 2 * X_mat[2, ])
  var_diff <- t(x_vec) %*% beta_vcov[var_names, var_names, drop = FALSE] %*% x_vec
  se_diff  <- sqrt(var_diff[[1]])
  se_diff  <- round(se_diff, 2)

  p_val <- round(2 * pnorm(abs(diff) / se_diff, lower.tail = FALSE), 2)

  res_df <- rbind(res_df,
                  data.frame(
                    slope1   = slope1,
                    slope2   = slope2,
                    diff     = diff,
                    se_diff  = se_diff,
                    p_value  = p_val
                  )
  )

  rownames(res_df) <- NULL
  return(res_df)
}

##' Wrapper for the baseline trajectory (beta_mu) -------------------------------
test_beta_mu_slope <- function(nknots, fit) {
  k_vec <- if (is.list(fit)) fit[[1]]$spline_inner_knots else fit$spline_inner_knots
  new_times <- c(k_vec[1] - 0.1, k_vec[1], k_vec[1] + 0.1)
  test_beta_slope(new_times, fit, var_name = "X0")
}

##' Wrapper for the treatment effect trajectory (beta_A) -----------------------
test_beta_A_slope <- function(nknots, fit) {
  k_vec <- if (is.list(fit)) fit[[1]]$spline_inner_knots else fit$spline_inner_knots
  new_times <- c(k_vec[1] - 0.1, k_vec[1], k_vec[1] + 0.1)
  test_beta_slope(new_times, fit, var_name = "A")
}
