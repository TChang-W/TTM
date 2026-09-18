##' Evaluate the estimated treatment effect trajectory on the retrospective time scale
##'
##' @description
##' Computes the estimated treatment effect trajectory,
##' \eqn{\beta_A(t^*)}, at user-supplied retrospective time points from a fitted
##' `TTM()` or `RetroJM()` object.
##'
##' The trajectory is formed from the treatment main effect together with the
##' treatment-by-spline interaction terms in the fitted retrospective-time
##' longitudinal model. Pointwise confidence intervals are constructed using a
##' normal approximation based on the estimated covariance matrix of the fixed
##' effects.
##'
##' @param new_times Numeric vector of retrospective time points
##'   \eqn{t^*} at which to evaluate the treatment effect trajectory.
##' @param results A fitted object returned by `TTM()` or `RetroJM()`.
##'   The object must contain spline knot information and the estimated fixed
##'   effects with their covariance matrix.
##' @param conf.level Numeric scalar giving the confidence level for pointwise
##'   confidence intervals. Must be strictly between `0` and `1`.
##'
##' @details
##' For a fitted retrospective-time model of the form
##'
##' \deqn{
##' E[Y(t^*)] = \beta_X^\top X +
##' \beta_{\mu}(t^*) + A \times \beta_A(t^*),
##' }
##'
##' this function evaluates \eqn{\beta_A(t^*)} at specified \eqn{t^*}
##'
##'
##' The spline basis is reconstructed using the boundary knots and internal knots
##' stored in `results`, so the supplied `new_times` should be on the same
##' retrospective-time scale as used when fitting the model.
##'
##' Confidence intervals are pointwise, not simultaneous.
##'
##' @return A named list with components:
##' \describe{
##'   \item{time_grid}{The input retrospective time points in `new_times`.}
##'   \item{estimate}{Estimated treatment effect values
##'   \eqn{\beta_A(t^*)}.}
##'   \item{standard_error}{Pointwise standard errors of the estimated treatment
##'   effect trajectory.}
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
##' trt_eff <- get_beta_A(
##'   new_times = seq(0, 10, by = 1),
##'   results = fit,
##'   conf.level = 0.95
##' )
##'
##' trt_eff$estimate
##' trt_eff$CI.lower
##' trt_eff$CI.upper
##' }
##'
##' @export
get_beta_A <- function(new_times, results, conf.level = 0.95) {
  res <- NULL
  if (conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be strictly between 0 and 1.")
  }
  if(results$method == "TTM"){

    if (is.null(results$weighted_gee)) {
      stop("`results` does not contain `weighted_gee`.")
    }

    Boundary.knots <- results$spline_boundary_knots
    knots <- results$spline_inner_knots

    beta_hat <- results$weighted_gee$fixed_effects
    beta_vcov <- results$weighted_gee$fixed_effects_vcov

    beta_A_vars <- c("A", grep("^AX_spline", names(beta_hat), value = TRUE))
    beta_A_coef <- beta_hat[beta_A_vars]

    if(results$type == "linear"){
      spline_val <- switch(results$spline_type,
                           "ns" = splines::ns(new_times, knots = knots,
                                              degree = 1, intercept = FALSE),
                           "bs" = splines::bs(new_times, knots = knots,
                                              degree = 1, intercept = FALSE)
      )}else{
        spline_val <- splines::ns(
          new_times,
          knots = knots,
          Boundary.knots = Boundary.knots,
          intercept = FALSE
        )}

    X_A <- cbind(A = 1, spline_val)
    beta_A <- as.vector(X_A %*% beta_A_coef)

    cov_beta_A <- beta_vcov[beta_A_vars, beta_A_vars, drop = FALSE]
    var_beta_A <- apply(X_A, 1, function(x) {
      as.numeric(t(x) %*% cov_beta_A %*% x)
    })
    sd_beta_A <- sqrt(var_beta_A)

    z_value <- -stats::qnorm((1 - conf.level) / 2)
    beta_A_lower <- beta_A - z_value * sd_beta_A
    beta_A_upper <- beta_A + z_value * sd_beta_A

    res <- list(
      time_grid = new_times,
      estimate = beta_A,
      standard_error = sd_beta_A,
      CI.lower = beta_A_lower,
      CI.upper = beta_A_upper
    )
  }
  if(results$method == "RetroJM"){
    name_map <- results$variable_name_map
    Boundary.knots <- results$spline_boundary_knots
    knots <- results$spline_inner_knots

    beta_lmer <- results$pooled_longitudinal_model$fixed_effects
    beta_lmer_vcov <- results$pooled_longitudinal_model$fixed_effects_vcov

    A_name <- name_map["A", "original_name"]
    beta_A_vars <- c(A_name, grep("^AX_spline", names(beta_lmer), value = TRUE))
    beta_A_coef <- beta_lmer[beta_A_vars]

    spline_val <- splines::ns(
      new_times,
      knots = knots,
      Boundary.knots = Boundary.knots
    )

    X_A <- cbind(`A` = 1, spline_val)
    beta_A <- as.vector(X_A %*% beta_A_coef)

    cov_beta_A <- beta_lmer_vcov[beta_A_vars, beta_A_vars, drop = FALSE]
    var_beta_A <- apply(X_A, 1, function(x) {
      as.numeric(t(x) %*% cov_beta_A %*% x)
    })
    sd_beta_A <- sqrt(var_beta_A)

    z_value <- -stats::qnorm((1 - conf.level) / 2)
    beta_A_lower <- beta_A - z_value * sd_beta_A
    beta_A_upper <- beta_A + z_value * sd_beta_A

    res <- list(
      time_grid = new_times,
      estimate = beta_A,
      standard_error = sd_beta_A,
      CI.lower = beta_A_lower,
      CI.upper = beta_A_upper
    )
  }
  return(res)
}
