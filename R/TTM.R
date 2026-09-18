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
##' @param outcome Character scalar giving the longitudinal outcome column name.
##' @param covariates Character vector of baseline covariate column names.
##' @param os_time Character scalar giving the observed survival time column
##'   name.
##' @param event Character scalar giving the event indicator column name.
##'   The function assumes `1 = death`, `2 = dropout`, and all other values are
##'   treated as neither death nor dropout.
##' @param t_pros Character scalar giving the prospective follow-up time column
##'   name.
##' @param jk_block_size Positive integer giving the number of delete-group
##'   jackknife groups used to estimate coefficient variability. Must be greater
##'   than 1.
##' @param n_inner_knot_list Positive integer vector giving the candidates of
##'   number of internal knots for the natural spline basis in retrospective
##'   time. The final value will be selected using criteria.
##' @param corstr Character scalar specifying the working correlation structure
##'   passed to `geepack::geeglm()`, such as `"independence"`, `"exchangeable"`,
##'   or `"ar1"`.
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
##' fit <- TTM(
##'   data_long = TTM_data,
##'   id = "id",
##'   treatment = "A",
##'   outcome = "Y",
##'   covariates = c("X1", "X2"),
##'   os_time = "OS_time",
##'   event = "event",
##'   t_pros = "t_pros",
##'   jk_block_size = NULL,
##'   n_inner_knot_list = range(10),
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
TTM <- function(data_long,
                id = "id",
                treatment = "A",
                outcome = "Y",
                covariates = c("X1", "X2"),
                os_time = "OS_time",
                event = "event",
                t_pros = "t_pros",
                jk_block_size = NULL,
                n_inner_knot_list = range(10),
                corstr = "ar1",
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


  if (!is.numeric(n_inner_knot_list) || any(n_inner_knot_list < 1) ) {
    stop("`n_inner_knot_list` must be a integer vector greater than or equal to 1.")
  }
  if(is.null(treatment)){
    treatment = "A"
    data_input$A = 0
  }

  data_input <- data_long
  data_input$death <- as.integer(data_input[[event]] == 1)
  data_input$dropout <- as.integer(data_input[[event]] == 2)

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
  knots_results <- list()

  for (i in seq_along(n_inner_knot_list)) {
    result <- try(
      TTM_knots(
        data_trial,
        n_inner_knot_list[i],
        corstr,
        jackknife = FALSE,
        jk_block_size = jk_block_size
      ),
      silent = TRUE
    )

    if (inherits(result, "try-error")) {
      message("Dropping after n_inner_knot = ", i, ": too many knots")
      break
    }

    knots_results[[i]] <- result
  }
  if(length(knots_results) > 0){
    keep <- sapply(knots_results, function(x) {
      !is.null(x$gee_model) &&
        !is.null(x$gee_model$converged) &&
        isTRUE(x$gee_model$converged)
    })

    knots_results <- knots_results[keep]

    gee_metrics <- gee_criteria(knots_results)
    n_inner_knot <- gee_metrics$n_inner_knot[which.min(gee_metrics[[criteria]])]

    final_res <- TTM_knots(data_trial, n_inner_knot, corstr,
                           jackknife = TRUE, jk_block_size = jk_block_size)
  }else{
    final_res <- TTM_knots(data_trial, n_inner_knot = 0, corstr,
                           jackknife = TRUE, jk_block_size = jk_block_size)
  }

  list(
    dropout_formula = final_res$dropout_formula,
    dropout_model = final_res$dropout_model,
    longitudinal_formula = final_res$longitudinal_formula,
    weighted_gee = final_res$gee_estimate,
    spline_boundary_knots = final_res$spline_boundary_knots,
    spline_inner_knots = final_res$spline_inner_knots,
    n_inner_knots = n_inner_knot,
    variable_name_map = name_map,
    weights = final_res$weights_df,
    gee_metrics = gee_metrics,
    data = data_trial,
    method = "TTM",
    type = "default"
  )
}

