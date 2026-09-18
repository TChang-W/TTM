##' Fit a weighted GEE terminal trend model on the retrospective time scale
##'
##' @description
##' Fits a retrospective-time longitudinal model using weighted generalized
##' estimating equations (WGEE).
##'
##' Retrospective time is the time from prospective follow-up time to the death
##' time. A natural spline basis is constructed for
##' retrospective time to allow flexible terminal-time trajectories. To account
##' for informative dropout and administrative censoring, subject-specific
##' weights are estimated. The weighted longitudinal mean model is
##' then fit with a Gaussian GEE among subjects with observed death events.
##'
##' Uncertainty for the regression coefficients is estimated using a
##' delete-group jackknife.
##'
##' @param data_long A long-format `data.frame` with one row per subject-visit.
##'   The data must include the columns specified by `id`, `treatment`,
##'   `outcome`, `covariates`, `os_time`, `event`, and `t_pros`, as well as
##'   `censoring_time` and `death_time`.
##' @param id Character scalar giving the subject identifier column name.
##' @param treatment Character scalar giving the treatment column name.
##'   Can only accept NULL or binary group indicator.
##' @param outcome Character scalar giving the longitudinal outcome column name.
##' @param covariates Character vector of baseline covariate column names.
##' @param os_time Character scalar giving the observed survival time column
##'   name.
##' @param event Character scalar giving the event indicator column name.
##'   The function assumes `1 = death`, `2 = dropout`, and all other values are
##'   treated as neither death nor dropout.
##' @param t_pros Character scalar giving the prospective follow-up time column
##'   name.
##' @param knots pre-specified knots. If null, the knot is selected automatically
##' @param jk_block_size Positive integer giving the number of delete-group
##'   jackknife groups used to estimate coefficient variability. Must be greater
##'   than 1.
##' @param corstr Character scalar specifying the working correlation structure
##'   passed to `geepack::geeglm()`, such as `"independence"`, `"exchangeable"`,
##'   or `"ar1"`.
##' @param spline_type Character scalar specifying the type of spline basis to
##'   use. Options are `"ns"` for natural splines or `"bs"` for B-splines.
##'   The default is `"bs"`.
##'
##' @details
##' The fitted longitudinal model includes:
##'
##' \deqn{
##' E[Y(t^*)] = \beta_{\mu}(t^*) + A \times \beta_A(t^*) + \beta_X^T X
##' }
##'
##' where \eqn{\beta_{\mu}(t^*)} and \eqn{\beta_A(t^*)} are represented by a
##' natural spline basis, and \eqn{t^* = D_i - t_{pros}} where \eqn{D_i} is the
##' death time.
##'
##' Dropout weights are derived from a Cox model for the hazard of dropout and
##' empirical estimates of censoring. The GEE fit is obtained with weights.
##' Only observations from subjects with `event == 1` are used in the final
##' longitudinal model fit.
##'
##' @return A named list with components:
##' \describe{
##'   \item{dropout_formula}{The formula used to fit the dropout Cox model.}
##'   \item{dropout_model}{A fitted `survival::coxph` object for dropout.}
##'   \item{longitudinal_formula}{The formula used in the weighted GEE fit.}
##'   \item{weighted_gee}{A list containing:
##'     \describe{
##'       \item{fixed_effects}{Jackknife mean of the GEE coefficient estimates.}
##'       \item{fixed_effects_vcov}{Delete-group jackknife covariance matrix for
##'       `fixed_effects`.}
##'       \item{fixed_effects_sd}{Standard errors from
##'       `fixed_effects_vcov`.}
##'       \item{jackknife_coefficients}{Matrix of leave-group-out coefficient
##'       estimates, one row per jackknife replicate.}
##'     }
##'   }
##'   \item{spline_boundary_knots}{Boundary knots used for the retrospective-time
##'   natural spline basis.}
##'   \item{spline_inner_knots}{Internal knots used for the retrospective-time
##'   natural spline basis.}
##'   \item{n_inner_knots}{The requested number of internal spline knots.}
##'   \item{variable_name_map}{Mapping between original variable names and the
##'   standardized names created during preprocessing.}
##'   \item{weights}{A `data.frame` containing subject-level analysis weights.}
##'   \item{method}{Character string identifying the fitted method, equal to
##'   `"TTM"`.}
##' }
##'
##'
##' @importFrom magrittr %>%
##' @importFrom survival Surv coxph basehaz
##' @importFrom splines ns
##' @importFrom geer geewa
##'
##' @examples
##' \dontrun{
##' library(TTM)
##' data("TTM_data", package = "TTM")
##'
##' fit <- TTM_linear(
##'   data_long = TTM_data,
##'   id = "id",
##'   treatment = NULL,
##'   outcome = "Y",
##'   covariates = c("X1", "X2"),
##'   os_time = "OS_time",
##'   event = "event",
##'   t_pros = "t_pros",
##'   knots = 6,
##'   jk_block_size = NULL,
##'   spline_type = "bs",
##'   corstr = "ar1"
##' )
##'
##' fit$weighted_gee$fixed_effects
##' fit$weighted_gee$fixed_effects_vcov
##'
##' plot_beta_mu(
##'   fit,
##'   xlim = c(0, stats::quantile(TTM_data$OS_time, 0.7)),
##'   conf.level = 0.95
##' )
##'
##' plot_beta_A(
##'   fit,
##'   xlim = c(0, stats::quantile(TTM_data$OS_time, 0.7)),
##'   conf.level = 0.95
##' )
##'
##' get_beta_mu(new_times = 1:10, TTM_results = fit, conf.level = 0.95)
##' get_beta_A(new_times = 1:10, TTM_results = fit, conf.level = 0.95)
##' }
##'
##' @export
TTM_linear <- function(data_long,
                       id = "id",
                       treatment = NULL,
                       outcome = "Y",
                       covariates = c("X1", "X2"),
                       os_time = "OS_time",
                       event = "event",
                       t_pros = "t_pros",
                       knots = NULL,
                       jk_block_size = NULL,
                       corstr = "ar1",
                       spline_type = "bs",
                       criteria = "RJC") {

  if (!is.data.frame(data_long)) {
    stop("`data_long` must be a data.frame.")
  }

  required_vars <- unique(c(
    id, treatment, outcome, covariates, os_time, event, t_pros
  ))
  missing_vars <- setdiff(required_vars, names(data_long))
  if (length(missing_vars) > 0) {
    stop("Missing required columns in `data_long`: ",
         paste(missing_vars, collapse = ", "))
  }


  data_input <- data_long
  data_input$death <- as.integer(data_input[[event]] == 1)
  data_input$dropout <- as.integer(data_input[[event]] == 2)
  if(is.null(treatment)){
    treatment = "A"
    data_input$A = 0
  }

  data_out <- clean_data(
    data_input = data_input,
    id = id,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    t_pros = t_pros,
    os_time = os_time,
    death = "death",
    dropout = "dropout"
  )

  data_trial <- data_out$data
  name_map <- data_out$name_map

  data_trial <- data_trial %>%
    dplyr::mutate(
      X0 = 1,
      t_retro = OS_time - t_pros,
      W1 = t_retro
    )

  final_res <- TTM_knots_linear(data_trial, corstr, knots, spline_type,
                                jackknife = TRUE, jk_block_size = jk_block_size)

  list(
    dropout_formula = final_res$dropout_formula,
    dropout_model = final_res$dropout_model,
    longitudinal_formula = final_res$longitudinal_formula,
    weighted_gee = final_res$gee_estimate,
    spline_boundary_knots = final_res$spline_boundary_knots,
    spline_inner_knots = final_res$spline_inner_knots,
    knots = knots,
    variable_name_map = name_map,
    weights = final_res$weights_df,
    data = data_trial,
    method = "TTM",
    type = "linear",
    spline_type = spline_type
  )
}

