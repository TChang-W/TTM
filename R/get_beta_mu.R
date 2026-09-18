##' Get the estimated baseline retrospective trajectory with pointwise confidence intervals
##'
##' Evaluates the estimated baseline longitudinal trajectory on the retrospective
##' timescale, denoted by \eqn{\beta_{\mu}(t^*)}, at user-supplied time points.
##' Pointwise normal approximation confidence intervals are returned based on
##' the weighted GEE fixed effect estimates from a fitted `TTM()` or `RetroJM()` object.
##'
##' The baseline trajectory corresponds to the mean spline-based effect of
##' retrospective time for the reference treatment group and baseline covariate
##' pattern encoded by the model intercept-like term `X0` and spline terms.
##'
##' @param new_times A numeric vector of retrospective time points at which the
##'   baseline trajectory is to be evaluated.
##' @param results A fitted object returned by `TTM()` or `RetroJM()`.
##' @param conf.level A numeric confidence level for the pointwise confidence
##'   interval. Must be strictly between `0` and `1`. The default is `0.95`.
##' @details
##' For a fitted retrospective-time model of the form
##'
##' \deqn{
##' E[Y(t^*)] = \beta_X^\top X +
##' \beta_{\mu}(t^*) + A \times \beta_A(t^*),
##' }
##'
##' this function evaluates \eqn{\beta_{\mu}(t^*)} at specified \eqn{t^*}
##'
##'
##' The spline basis is reconstructed using the boundary knots and internal knots
##' stored in `results`, so the supplied `new_times` should be on the same
##' retrospective-time scale as used when fitting the model.
##'
##' Confidence intervals are pointwise, not simultaneous.
##' @return A list containing:
##' \describe{
##'   \item{time_grid}{The user-supplied retrospective time points.}
##'   \item{estimate}{The estimated baseline trajectory \eqn{\beta_{\mu}(t^*)}.}
##'   \item{standard_error}{The pointwise standard error of the estimated trajectory.}
##'   \item{CI.lower}{The lower bound of the pointwise confidence interval.}
##'   \item{CI.upper}{The upper bound of the pointwise confidence interval.}
##' }
##'
##' @export
get_beta_mu <- function(new_times, results, conf.level = 0.95) {
  res <- NULL
  if(results$method == "TTM"){
    if (conf.level <= 0 || conf.level >= 1) {
      stop("`conf.level` must be strictly between 0 and 1.")
    }

    if (is.null(results$weighted_gee)) {
      stop("`results` does not contain `weighted_gee`.")
    }

    Boundary.knots <- results$spline_boundary_knots
    knots <- results$spline_inner_knots

    beta_hat <- results$weighted_gee$fixed_effects
    beta_vcov <- results$weighted_gee$fixed_effects_vcov

    beta_mu_vars <- c("X0", grep("^X_spline", names(beta_hat), value = TRUE))
    beta_mu_coef <- beta_hat[beta_mu_vars]

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

    res <- list(
      time_grid = new_times,
      estimate = beta_mu,
      standard_error = sd_beta_mu,
      CI.lower = beta_mu_lower,
      CI.upper = beta_mu_upper
    )
  }
  if(results$method == "RetroJM"){
    Boundary.knots <- results$spline_boundary_knots
    knots <- results$spline_inner_knots

    beta_lmer <- results$pooled_longitudinal_model$fixed_effects
    beta_lmer_vcov <- results$pooled_longitudinal_model$fixed_effects_vcov

    beta_mu_vars <- c("(Intercept)", grep("^X_spline", names(beta_lmer), value = TRUE))
    beta_mu_coef <- beta_lmer[beta_mu_vars]

    spline_val <- splines::ns(
      new_times,
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

    res <- list(
      time_grid = new_times,
      estimate = beta_mu,
      standard_error = sd_beta_mu,
      CI.lower = beta_mu_lower,
      CI.upper = beta_mu_upper
    )
  }
  return(res)
}
