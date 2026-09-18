##' Fit the two-stage model for survival imputation and longitudinal analysis
##'
##' @description
##' Fits a two-stage retrospective joint modeling procedure. In stage 1, a Cox
##' proportional hazards model is fit to overall survival and used to impute death
##' times for censored individuals. In stage 2, the imputed datasets are used to
##' fit a linear mixed effects model for the longitudinal outcome measured on a
##' retrospective timescale from death. Fixed effects estimates and associated
##' uncertainty are then pooled across imputations using standard multiple
##' imputation combining rules.
##'
##' @param data_long A long-format `data.frame` containing one row per subject per visit.
##'   The dataset must include the subject identifier, treatment variable,
##'   longitudinal outcome, survival time-to-event outcome, prospective follow-up time, and any
##'   baseline covariates used in the models.
##'
##' @param id A character string giving the name of the subject identifier column.
##'   Each unique value should correspond to one individual.
##'
##' @param treatment A character string giving the name of the treatment or exposure
##'   variable. This variable is included in both the survival and longitudinal
##'   models. It is typically binary, but other codings may be used if supported by
##'   the model specification.
##'
##' @param outcome A character string giving the name of the longitudinal outcome
##'   variable to be modeled, such as a quality-of-life score.
##'
##' @param covariates A character vector of baseline covariate names to include as
##'   adjustment variables in both the survival and longitudinal submodels.
##'
##' @param os_time A character string giving the name of the overall survival time
##'   variable. This should represent observed event or censoring time on the
##'   prospective timescale.
##'
##' @param death A character string giving the name of the event indicator variable.
##'   This should typically be coded as `1` for death and `0` for censoring.
##'
##' @param t_pros A character string giving the name of the prospective follow-up
##'   time variable for each longitudinal observation. This is used together with
##'   survival time to define retrospective time from death or imputed death.
##'
##' @param random_effects A character vector specifying the random effects terms to
##'   include in the linear mixed effects model. For example, `c("1", "A")`
##'   specifies a random intercept and a random slope for treatment.
##'
##' @param n_impute A positive integer giving the number of imputed death times, and
##'   therefore the number of imputed datasets, to generate for censored subjects.
##'   Larger values generally improve stability of pooled estimates at the cost of
##'   increased computation time.
##'
##' @param n_inner_knot A positive integer giving the number of internal knots used
##'   for the natural spline basis that models the retrospective time trend in the
##'   longitudinal submodel. Increased number of knots may improve the model capacity
##'   to capture complicated trajectories at the risk of capturing noises.
##'
##' @details
##'
##' This function implements a two-stage approach for analyzing longitudinal
##' outcomes that are conceptually aligned by time before death.
##'
##' Cox proportional hazards model:
##' \eqn{\lambda_i(t)=\exp{\Big(\alpha_0(t)+A_i\alpha_A(t)\Big)}}
##'
##' Longitudinal linear mixed effect Model:
##' \eqn{Y_i(t^*)=\beta_{\mu}(t^*)+A_i\beta_A(t^*)+X_i^T\psi_X+W_ib_i+\epsilon_i(t^*)}
##'
##' In the first stage, a Cox proportional hazards model is fit using treatment and
##' baseline covariates. For subjects who are censored, conditional survival
##' distributions beyond the censoring time are derived from the fitted Cox model,
##' and death times are multiply imputed.
##'
##' In the second stage, each imputed dataset is reconstructed by replacing
##' censored survival times with imputed death times and fit the longitudinal
##' process on the retrospective timescale. Then a natural spline basis
##' is constructed on the retrospective time scale to capture the baseline
##' trajectory and treatment caused trajectory for outcome  variable, and a
##' linear mixed effects model is fit to the longitudinal outcome.
##'
##' The longitudinal model includes treatment, baseline covariates, spline terms
##' for retrospective time, and treatment-by-spline interactions. Random effects
##' are included according to `random_effects` by user specification.
##'
##' Fixed effects estimates and covariance matrices from the mixed models are pooled
##' across imputations using Rubin-style combining rules. Variance components are
##' summarized by simple averaging across imputations.
##'
##' This function assumes that:
##' \itemize{
##'   \item `data_long` is in long format, with repeated outcome measurements per subject,
##'   \item survival and censoring information are available at the subject level,
##'   \item the longitudinal outcome is observed prior to death or censoring,
##'   \item covariates included in the models are measured consistently across rows.
##' }
##'
##' @return A named list with components:
##' \describe{
##'   \item{survival_formula}{The survival model formula used in the Cox model.}
##'   \item{survival_model}{The fitted `survival::coxph` object from stage 1.}
##'   \item{longitudinal_formula}{The linear mixed effects model formula used in stage 2.}
##'   \item{imputed_death_times}{A matrix of imputed death times for censored subjects,
##'   with one column per imputation and one row per censored subject.}
##'   \item{pooled_longitudinal_model}{A list containing pooled results from the
##'   longitudinal mixed effects models:
##'   \describe{
##'     \item{fixed_effects}{Pooled fixed effects estimates.}
##'     \item{fixed_effects_vcov}{Pooled covariance matrix of the fixed effects estimates.}
##'     \item{residual_variance}{Average residual variance across imputations.}
##'     \item{random_effects_covariance}{Average random effects covariance parameters
##'     across imputations.}
##'   }}
##'   \item{spline_boundary_knots}{Boundary knots used for the natural spline basis.}
##'   \item{spline_inner_knots}{Internal knots used for the natural spline basis.}
##'   \item{n_inner_knots}{The number of internal knots used.}
##'   \item{variable_name_map}{A lookup table showing original variable names and any
##'   internally transformed names used during preprocessing.}
##' }
##'
##' @importFrom magrittr %>%
##' @importFrom ggplot2 aes
##' @importFrom utils head tail
##' @importFrom survival Surv survSplit coxph basehaz
##'
##' @examples
##' \dontrun{
##' library(TTM)
##' data("Retro_data", package = "TTM")
##'
##' fit <- TTM(
##'   data_long = Retro_data,
##'   id = "id",
##'   treatment = "A",
##'   outcome = "Y",
##'   covariates = c("X1", "X2"),
##'   os_time = "OS_time",
##'   death = "death",
##'   t_pros = "t_pros",
##'   random_effects = c("1", "A"),
##'   n_impute = 10,
##'   n_inner_knot = 3
##' )
##'
##' fit$pooled_longitudinal_model$fixed_effects
##' fit$pooled_longitudinal_model$fixed_effects_vcov
##'
##' ##' plot_beta_mu(fit, xlim = c(0, quantile(Retro_data$OS_time, 0.7)), conf.level = 0.95)
##' plot_beta_A(fit, xlim = c(0, quantile(Retro_data$OS_time, 0.7)), conf.level = 0.95)
##'
##' get_beta_mu(new_times = 1:10, TTM_results = fit, conf.level = 0.95)
##' get_beta_A(new_times = 1:10, TTM_results = fit, conf.level = 0.95)
##' }
##'
##'
##' @export
RetroJM <- function(data_long,
                    id = "id",
                    treatment = "A",
                    outcome = "Y",
                    covariates = c("X1", "X2"),
                    os_time = "OS_time",
                    death = "death",
                    t_pros = "t_pros",
                    random_effects = c("1", "A"),
                    jackknife_size = 50,
                    n_inner_knot = 3) {
  data_out = clean_data(data_long,
                        id = id,
                        outcome = outcome,
                        treatment = treatment,
                        covariates = covariates,
                        os_time = os_time,
                        death = death)
  data_cleaned = data_out$data
  name_map = data_out$name_map
  data_cleaned <- data_cleaned %>%
    dplyr::mutate(
      t_retro = OS_time - t_pros
    )
  time_grid <- data_cleaned$t_retro
  knots <- stats::quantile(time_grid, (1:n_inner_knot) / (n_inner_knot + 1))
  Boundary.knots <- range(time_grid)
  data_spline <- splines::ns(
    time_grid,
    knots = knots,
    Boundary.knots = Boundary.knots,
    intercept = FALSE
  )
  X_spline <- data_spline
  colnames(X_spline) <- paste0("X_spline", seq_len(ncol(data_spline)))
  X_spline_vars <- colnames(X_spline)
  data_cleaned <- cbind(data_cleaned, X_spline)
  for (xv in X_spline_vars) {
    data_cleaned[[paste0("A", xv)]] <- data_cleaned$A * data_cleaned[[xv]]
  }
  data_combined <- data_cleaned
  data_death <- dplyr::filter(data_cleaned, death == 1)
  data_censor <- dplyr::filter(data_cleaned, death != 1)
  surv_vars = name_map[c(id, treatment, covariates, os_time, death), "transformed_name"]
  random_effects_mapped <- name_map[random_effects, "transformed_name"]
  random_effects_mapped <- ifelse(is.na(random_effects_mapped), random_effects, random_effects_mapped)
  data_predict <- unique(
    data_censor[, surv_vars]
  )
  Surv <- survival::Surv
  surv.formula <- stats::as.formula(
    paste0("Surv(", os_time, ", ", death, ") ~ 1 + ",
           treatment, " + ", paste0(covariates, collapse = " + "))
  )
  data_combined_cut <- survival::survSplit(
    surv.formula,
    data_combined,
    cut = unique(data_combined$OS_time[data_combined$death == 1]),
    episode = "t_cut"
  ) %>%
    unique()
  data_combined_cut <- data_combined_cut[
    (data_combined_cut$OS_time - data_combined_cut$tstart > 1e-03),
  ]
  surv_model <- survival::coxph(
    surv.formula,
    data = data_combined_cut
  )
  base_cumhaz_df <- survival::basehaz(surv_model, centered = FALSE)
  base_cumhaz_unique <- base_cumhaz_df[!duplicated(base_cumhaz_df$hazard), ]
  base_insthaz_df <- base_cumhaz_unique
  base_insthaz_df$hazard <- diff(c(0, base_cumhaz_unique$hazard)) /
    diff(c(0, base_cumhaz_unique$time))
  base_insthaz_df$time <- c(0, head(base_insthaz_df$time, -1))
  death_imputed <- matrix(NA_real_, nrow = nrow(data_predict), ncol = n_impute)
  risk <- rep(NA_real_, nrow(data_predict))
  for (i in seq_len(nrow(data_predict))) {
    risk[i] <- stats::predict(
      surv_model,
      newdata = data_predict[i, ],
      type = "risk"
    )
    pt_insthaz_df <- base_insthaz_df %>%
      dplyr::mutate(
        hazard = hazard * risk[i]
      )
    pt_cumhaz_df <- data.frame(
      time = pt_insthaz_df$time,
      hazard = cumsum(
        head(c(0, pt_insthaz_df$hazard), -1) * c(0, diff(pt_insthaz_df$time))
      )
    )
    tc <- data_predict$OS_time[i]
    cumhaz_tc_i <- stats::approx(
      x = pt_cumhaz_df$time,
      y = pt_cumhaz_df$hazard,
      xout = tc,
      method = "linear",
      rule = 2
    )$y
    cumhaz_i <- if (is.na(cumhaz_tc_i)) 0 else pmax(pt_cumhaz_df$hazard - cumhaz_tc_i, 0)
    surv_prob_i <- stats::runif(n_impute)
    for (j in seq_len(n_impute)) {
      target <- -log(surv_prob_i[j])
      cumhaz_unique_idx <- !duplicated(cumhaz_i)
      cumhaz_unique <- cumhaz_i[cumhaz_unique_idx]
      time_unique <- pt_cumhaz_df$time[cumhaz_unique_idx]

      if (target <= max(cumhaz_unique)) {
        death_imputed[i, j] <- stats::approx(
          x = cumhaz_unique,
          y = time_unique,
          xout = target,
          method = "linear",
          rule = 2
        )$y
      } else {
        death_imputed[i, j] <- max(time_unique) +
          (target - max(cumhaz_unique)) / tail(pt_insthaz_df$hazard, 1)
      }
    }
  }
  beta_lmer_impute <- NULL
  tau2_lmer_impute <- NULL
  Sigma_lmer_impute <- NULL
  beta_lmer_vcov_impute <- list()
  for (i in seq_len(n_impute)) {
    ED <- data.frame(id = data_predict$id, death_imputed = death_imputed[, i])
    ED <- dplyr::left_join(data_censor, ED, by = "id") %>%
      dplyr::mutate(
        death_time = death_imputed,
        OS_time = death_imputed,
        death = 1
      ) %>%
      dplyr::select(-death_imputed)
    data_combined <- rbind(data_death, ED) %>%
      dplyr::mutate(
        t_retro_imputed = OS_time - t_pros,
        OS_time_sqrt = sqrt(OS_time)
      )
    time_grid <- data_combined$t_retro_imputed
    data_combined_Boundary.knots <- range(death_imputed)
    data_combined_knots = quantile(death_imputed, (1:n_inner_knot) / (n_inner_knot + 1), na.rm = T)

    data_spline <- splines::ns(
      time_grid,
      knots = data_combined_knots,
      Boundary.knots = data_combined_Boundary.knots,
      intercept = FALSE
    )
    X_spline_imputed <- data_spline
    colnames(X_spline_imputed) <- paste0("X_spline", seq_len(ncol(X_spline_imputed)))
    X_spline_vars <- colnames(X_spline_imputed)
    data_combined[, X_spline_vars] <- X_spline_imputed
    for (xv in X_spline_vars) {
      data_combined[[paste0("A", xv)]] <- data_combined$A * data_combined[[xv]]
    }
    lmer_formula <- stats::as.formula(
      paste(
        "Y ~ ", treatment,
        "+", paste0(covariates, collapse = "+"), "+",
        paste(X_spline_vars, collapse = "+"),
        "+",
        paste(paste0("A", X_spline_vars), collapse = " + "),
        # "+(", paste0(random_effects_mapped, collapse = "+"), "| id)",
        "+ (1 + A||id)"
      )
    )
    lmer_model <- lme4::lmer(lmer_formula, REML = FALSE, data = data_combined)
    beta_lmer_impute <- rbind(beta_lmer_impute, lme4::fixef(lmer_model))
    beta_lmer_vcov_impute[[i]] <- stats::vcov(lmer_model)
    tau2_lmer_impute <- c(tau2_lmer_impute, attr(lme4::VarCorr(lmer_model), "sc")^2)
    Sigma_lmer_impute <- rbind(
      Sigma_lmer_impute,
      c(lapply(lme4::VarCorr(lmer_model), as.matrix)$id)
    )
  }
  beta_lmer <- colMeans(beta_lmer_impute)
  tau2_lmer <- mean(tau2_lmer_impute)
  Sigma_lmer <- colMeans(Sigma_lmer_impute)
  beta_U_lmer <- Reduce(`+`, beta_lmer_vcov_impute) / n_impute
  beta_diff_lmer <- sweep(beta_lmer_impute, 2, beta_lmer, "-")
  beta_B_lmer <- crossprod(beta_diff_lmer) / (n_impute - 1)
  beta_lmer_vcov <- beta_U_lmer + (1 + 1 / n_impute) * beta_B_lmer

  list(
    survival_formula = surv.formula,
    survival_model = surv_model,
    longitudinal_formula = lmer_formula,
    imputed_death_times = death_imputed,
    pooled_longitudinal_model = list(
      fixed_effects = beta_lmer,
      fixed_effects_vcov = beta_lmer_vcov,
      residual_variance = tau2_lmer,
      random_effects_covariance = Sigma_lmer
    ),
    spline_boundary_knots = data_combined_Boundary.knots,
    spline_inner_knots = data_combined_knots,
    n_inner_knots = n_inner_knot,
    variable_name_map = name_map,
    method = "RetroJM"
  )
}
