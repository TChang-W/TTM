##' Plot the estimated baseline trajectory on the retrospective time scale
##'
##' @description
##' Plots the estimated baseline retrospective trajectory,
##' \eqn{\beta_{\mu}(t^*)}, together with pointwise confidence intervals from a
##' fitted `TTM()` or `RetroJM()` object.
##'
##' The plotted trajectory is formed from the intercept-like baseline term
##' together with the spline terms in the fitted retrospective-time
##' longitudinal model. Pointwise confidence intervals are constructed using a
##' normal approximation based on the estimated covariance matrix of the fixed
##' effects.
##'
##' @param results A fitted object returned by `TTM()` or `RetroJM()`.
##'   The object must contain spline knot information and the estimated fixed
##'   effects with their covariance matrix.
##' @param xlim Optional numeric vector of length 2 giving the range of
##'   retrospective time to display. If `NULL`, the spline boundary range stored
##'   in `results` is used.
##' @param conf.level Numeric scalar giving the confidence level for pointwise
##'   confidence intervals. Must be strictly between `0` and `1`.
##'
##' @details
##' This function evaluates the baseline trajectory over a fine grid of
##' retrospective time values and plots the estimated curve together with its
##' lower and upper pointwise confidence limits.
##'
##' The spline basis is reconstructed using the boundary knots and internal knots
##' stored in `results`, so the plotted time axis should be interpreted on the
##' same retrospective-time scale used when fitting the model.
##'
##' Confidence intervals are pointwise, not simultaneous.
##'
##' @return Invisibly returns a named list with components:
##' \describe{
##'   \item{plot}{The `ggplot2` plot object.}
##'   \item{time_grid}{The retrospective time grid used for plotting.}
##'   \item{estimate}{Estimated baseline trajectory values
##'   \eqn{\beta_{\mu}(t^*)}.}
##'   \item{standard_error}{Pointwise standard errors of the estimated baseline
##'   trajectory.}
##'   \item{CI.lower}{Lower pointwise confidence limits.}
##'   \item{CI.upper}{Upper pointwise confidence limits.}
##' }
##'
##' @seealso
##' `TTM()`, `get_beta_mu()`, `plot_beta_A()`
##'
##' @examples
##' \dontrun{
##' fit <- TTM(
##'   data_long = TTM_data,
##'   id = "id",
##'   treatment = "A",
##'   outcome = "Y",
##'   covariates = c("X1", "X2"),
##'   os_time = "OS_time",
##'   event = "event",
##'   t_pros = "t_pros"
##' )
##'
##' plot_beta_mu(
##'   results = fit,
##'   xlim = c(0, 10),
##'   conf.level = 0.95
##' )
##' }
##' @export
plot_beta_mu <- function(results, xlim = NULL, conf.level = 0.95) {
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be strictly between 0 and 1.")
  }
  if(results$method == "TTM"){

    Boundary.knots <- results$spline_boundary_knots
    knots <- results$spline_inner_knots

    if (is.null(xlim)) {
      a <- Boundary.knots[1]
      b <- Boundary.knots[2]
    } else {
      if (!is.numeric(xlim) || length(xlim) != 2) {
        stop("xlim must be a numeric vector of length 2.")
      }
      a <- xlim[1]
      b <- xlim[2]
    }

    t_val <- seq(from = a, to = b, length.out = 1000)

    beta_hat <- results$weighted_gee$fixed_effects
    beta_vcov <- results$weighted_gee$fixed_effects_vcov

    beta_mu_vars <- c("X0", grep("^X_spline", names(beta_hat), value = TRUE))
    beta_mu_coef <- beta_hat[beta_mu_vars]

    if(results$type == "linear"){
      spline_val <- switch(results$spline_type,
                           "ns" = splines::ns(t_val, knots = knots,
                                              intercept = FALSE),
                           "bs" = splines::bs(t_val, knots = knots,
                                              degree = 1, intercept = FALSE)
      )}else{
        spline_val <- splines::ns(
          t_val,
          knots = knots,
          Boundary.knots = Boundary.knots,
          intercept = FALSE
        )}

    X_mu <- cbind(X0 = 1, spline_val)
    beta_mu <- as.vector(X_mu %*% beta_mu_coef)

    cov_beta_mu <- beta_vcov[beta_mu_vars, beta_mu_vars, drop = FALSE]
    var_beta_mu <- apply(X_mu, 1, function(x) {
      as.numeric(t(x) %*% cov_beta_mu %*% x)
    })
    sd_beta_mu <- sqrt(var_beta_mu)

    z_value <- -stats::qnorm((1 - conf.level) / 2)
    beta_mu_lower <- beta_mu - z_value * sd_beta_mu
    beta_mu_upper <- beta_mu + z_value * sd_beta_mu

    plot_df <- data.frame(
      time_grid = t_val,
      estimate = beta_mu,
      lower = beta_mu_lower,
      upper = beta_mu_upper
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time_grid, y = estimate)) +
      ggplot2::geom_line(color = "blue", linewidth = 1) +
      ggplot2::geom_line(
        ggplot2::aes(y = lower),
        color = "blue",
        linetype = "dashed"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = upper),
        color = "blue",
        linetype = "dashed"
      ) +
      ggplot2::labs(
        x = "Retrospective Time (t*)",
        y = expression(beta[mu](t^"*")),
        title = "Estimated Baseline Retrospective Trajectory with Pointwise Confidence Interval"
      ) +
      ggplot2::theme_minimal()

    print(p)

    return(invisible(list(
      plot = p,
      time_grid = t_val,
      estimate = beta_mu,
      standard_error = sd_beta_mu,
      CI.lower = beta_mu_lower,
      CI.upper = beta_mu_upper
    )))
  }
  if(results$method == "RetroJM"){
    Boundary.knots <- RetroJM_results$spline_boundary_knots
    knots <- RetroJM_results$spline_inner_knots
    if (is.null(xlim)) {
      a <- Boundary.knots[1]
      b <- tail(knots, 1)
    } else {
      if (!is.numeric(xlim) || length(xlim) != 2) {
        stop("xlim must be a numeric vector of length 2.")
      }
      a <- xlim[1]
      b <- xlim[2]
    }

    t_val <- seq(
      from = a,
      to = b,
      length.out = 1000
    )

    beta_lmer <- RetroJM_results$pooled_longitudinal_model$fixed_effects
    beta_lmer_vcov <- RetroJM_results$pooled_longitudinal_model$fixed_effects_vcov

    beta_mu_vars <- c("(Intercept)", grep("^X_spline", names(beta_lmer), value = TRUE))
    beta_mu_coef <- beta_lmer[beta_mu_vars]

    spline_val <- splines::ns(
      t_val,
      knots = knots,
      Boundary.knots = Boundary.knots
    )

    X_mu <- cbind(`(Intercept)` = 1, spline_val)
    beta_mu <- as.vector(X_mu %*% beta_mu_coef)

    cov_beta_mu <- beta_lmer_vcov[beta_mu_vars, beta_mu_vars, drop = FALSE]
    var_beta_mu <- apply(X_mu, 1, function(x) {
      as.numeric(t(x) %*% cov_beta_mu %*% x)
    })
    sd_beta_mu <- sqrt(var_beta_mu)

    z_value <- -stats::qnorm((1 - conf.level) / 2)
    beta_mu_lower <- beta_mu - z_value * sd_beta_mu
    beta_mu_upper <- beta_mu + z_value * sd_beta_mu

    plot_df <- data.frame(
      time_grid = t_val,
      estimate = beta_mu,
      lower = beta_mu_lower,
      upper = beta_mu_upper
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time_grid, y = estimate)) +
      ggplot2::geom_line(color = "blue", linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(y = lower), color = "blue", linetype = "dashed") +
      ggplot2::geom_line(ggplot2::aes(y = upper), color = "blue", linetype = "dashed") +
      ggplot2::labs(
        x = "Retrospective Time (t*)",
        y = expression(beta[mu](t^"*")),
        title = "Estimated Baseline Retrospective Trajectory with Pointwise Confidence Interval"
      ) +
      ggplot2::theme_minimal()

    print(p)

    invisible(list(
      plot = p,
      time_grid = t_val,
      estimate = beta_mu,
      standard_error = sd_beta_mu,
      CI.lower = beta_mu_lower,
      CI.upper = beta_mu_upper
    ))
  }
}
