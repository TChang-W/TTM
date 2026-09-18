TTM_knots <- function(dat, n_inner_knot, corstr,
                      jackknife = FALSE, jk_block_size){
  time_grid <- dat$t_retro
  if(n_inner_knot > 0){
  knots <- stats::quantile(
    time_grid,
    probs = (1:n_inner_knot) / (n_inner_knot + 1),
    na.rm = TRUE
  )
  Boundary.knots <- range(time_grid, na.rm = TRUE)

  data_spline <- splines::ns(
    time_grid,
    knots = knots,
    Boundary.knots = Boundary.knots,
    intercept = FALSE
  )

  X_spline <- data_spline
  colnames(X_spline) <- paste0("X_spline", seq_len(ncol(data_spline)))
  X_spline_vars <- colnames(X_spline)
  dat <- cbind(dat, X_spline)

  for (xv in X_spline_vars) {
    dat[[paste0("A", xv)]] <- dat$A * dat[[xv]]
  }

  covariate_terms <- grep("^X[0-9]+$", names(dat), value = TRUE)
  covariate_terms <- setdiff(covariate_terms, "X0")

  # test if A is valid group indicator
  if (length(unique(na.omit(dat$A))) > 1) {
    AX_spline_vars <- paste0("A", X_spline_vars)
    X_var <- unique(c("X0", "A", covariate_terms, X_spline_vars, paste0("A", X_spline_vars)))
    X_var <- X_var[X_var %in% names(dat)]

    dropout_formula <- stats::as.formula(
      paste0(
        "survival::Surv(OS_time, dropout) ~ 1 + A",
        if (length(covariate_terms) > 0) {
          paste0(" + ", paste(covariate_terms, collapse = " + "))
        } else {
          ""
        }
      )
    )
    longitudinal_formula <- stats::as.formula(
      paste(
        "Y ~ -1 + X0 + A",
        if (length(covariate_terms) > 0) paste("+", paste(covariate_terms, collapse = " + ")) else "",
        if (length(X_spline_vars) > 0) paste("+", paste(X_spline_vars, collapse = " + ")) else "",
        if (length(AX_spline_vars) > 0) paste("+", paste(AX_spline_vars, collapse = " + ")) else ""
      )
    )

  }else{
    AX_spline_vars <- NULL
    X_var <- unique(c("X0", covariate_terms, X_spline_vars))
    X_var <- X_var[X_var %in% names(dat)]
    dropout_formula <- stats::as.formula(
      paste0(
        "survival::Surv(OS_time, dropout) ~ 1",
        if (length(covariate_terms) > 0) {
          paste0(" + ", paste(covariate_terms, collapse = " + "))
        } else {
          ""
        }
      )
    )
    longitudinal_formula <- stats::as.formula(
      paste(
        "Y ~ -1 + X0",
        if (length(covariate_terms) > 0) paste("+", paste(covariate_terms, collapse = " + ")) else "",
        if (length(X_spline_vars) > 0) paste("+", paste(X_spline_vars, collapse = " + ")) else "",
        if (length(AX_spline_vars) > 0) paste("+", paste(AX_spline_vars, collapse = " + ")) else ""
      )
    )
  }



  }else{
    covariate_terms <- grep("^X[0-9]+$", names(dat), value = TRUE)
    covariate_terms <- setdiff(covariate_terms, "X0")
    if (length(unique(na.omit(dat$A))) > 1) {
      X_var <- unique(c("X0", "A", covariate_terms))
      X_var <- X_var[X_var %in% names(dat)]

      dropout_formula <- stats::as.formula(
        paste0(
          "survival::Surv(OS_time, dropout) ~ A",
          if (length(covariate_terms) > 0) {
            paste0(" + ", paste(covariate_terms, collapse = " + "))
          } else {
            ""
          }
        )
      )
      longitudinal_formula <- stats::as.formula(
        paste(
          "Y ~ -1 + X0 + A",
          if (length(covariate_terms) > 0) paste("+", paste(covariate_terms, collapse = " + ")) else ""
        )
      )
    }else{
      X_var <- unique(c("X0", covariate_terms))
      X_var <- X_var[X_var %in% names(dat)]

      dropout_formula <- stats::as.formula(
        paste0(
          "survival::Surv(OS_time, dropout) ~ 1",
          if (length(covariate_terms) > 0) {
            paste0(" + ", paste(covariate_terms, collapse = " + "))
          } else {
            ""
          }
        )
      )
      longitudinal_formula <- stats::as.formula(
        paste(
          "Y ~ -1 + X0",
          if (length(covariate_terms) > 0) paste("+", paste(covariate_terms, collapse = " + ")) else ""
        )
      )
    }

  }

  weights_df <- get_weights(dat, dropout_formula)
  dat <- dplyr::left_join(dat, weights_df, by = "id")
  data_death_w <- dat[dat$death == 1 & dat$weights > 0, , drop = FALSE]

  if (nrow(data_death_w) == 0) {
    stop("No death observations found after preprocessing.")
  }

  if (any(is.na(data_death_w$weights)) || any(data_death_w$weights <= 0)) {
    stop("Nonpositive or missing weights detected among death observations.")
  }

  scale_vars <- unique(c(X_var))
  scale_vars <- scale_vars[scale_vars %in% names(data_death_w)]

  data_death_w[, scale_vars] <- sqrt(1 / data_death_w$weights) * data_death_w[, scale_vars]
  data_death_w[, "Y"] <- sqrt(1 / data_death_w$weights) * data_death_w[, "Y"]


  gee_model <- geewa(
    formula = longitudinal_formula,
    family = gaussian(link = "identity"),
    data = data_death_w,
    id = id,
    corstr = corstr,
    method = "gee"
  )

  if(jackknife == TRUE){
    uniq_id <- unique(data_death_w$id)
  if(is.null(jk_block_size)){
    jk_size <- select_jk_block_size(length(uniq_id))
  }else{
    jk_size <- as.integer(jk_block_size)
  }

  id_jk <- jk_split(uniq_id, jk_size)

  beta_gee_weighted_jk <- NULL
  K <- length(id_jk)

  for (i in seq_len(K)) {
    leaveout_id <- unique(id_jk[[i]])
    data_death_w_jk <- data_death_w[!(data_death_w$id %in% leaveout_id), , drop = FALSE]

    if (nrow(data_death_w_jk) == 0) next

    fit <- try(
      geepack::geeglm(
        longitudinal_formula,
        family = gaussian,
        id = data_death_w_jk$id,
        data = data_death_w_jk,
        corstr = corstr
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) next

    beta_gee_weighted_jk <- rbind(beta_gee_weighted_jk, fit$coefficients)
  }

  if (is.null(beta_gee_weighted_jk) || nrow(beta_gee_weighted_jk) < 2) {
    stop("Jackknife failed to produce enough valid replicates.")
  }

  beta_gee_weighted_mean <- apply(beta_gee_weighted_jk, 2, mean)
  beta_gee_weighted_centered <- sweep(beta_gee_weighted_jk, 2, beta_gee_weighted_mean, "-")
  beta_gee_weighted_vcov <- (K - 1) / K * (t(beta_gee_weighted_centered) %*% beta_gee_weighted_centered)
  beta_gee_weighted_sd <- sqrt(diag(beta_gee_weighted_vcov))
  gee_estimate = list(
    fixed_effects = beta_gee_weighted_mean,
    fixed_effects_vcov = beta_gee_weighted_vcov,
    fixed_effects_sd = beta_gee_weighted_sd
  )
  }else{
    gee_estimate <- NULL
  }

  list(
    dropout_formula = dropout_formula,
    longitudinal_formula = longitudinal_formula,
    weights_df = weights_df,
    gee_model = gee_model,
    spline_boundary_knots = attr(data_spline, "Boundary.knots"),
    spline_inner_knots = attr(data_spline, "knots"),
    n_inner_knots = n_inner_knot,
    gee_estimate = gee_estimate
  )

}
