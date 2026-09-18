clean_data <- function(data_input,
                       id = "id",
                       treatment = "A",
                       outcome = "Y",
                       covariates = c("X1", "X2"),
                       t_pros = "t_pros",
                       os_time = "OS_time",
                       death = "death",
                       dropout = "dropout") {
  stopifnot(is.data.frame(data_input))

  needed <- c(id, treatment, outcome, covariates, t_pros, os_time, death, dropout)
  missing_cols <- setdiff(needed, names(data_input))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- data_input
  if(is.null(covariates)){
    mapped_names <- c(
      "id", "A", "Y",
      "t_pros", "OS_time", "death", "dropout"
    )
  }else{
    mapped_names <- c(
      "id", "A", "Y",
      paste0("X", seq_along(covariates)),
      "t_pros", "OS_time", "death", "dropout"
    )
  }


  original_names <- c(
    id, treatment, outcome,
    covariates, t_pros, os_time, death, dropout
  )

  name_map <- data.frame(
    original_name = original_names,
    transformed_name = mapped_names,
    stringsAsFactors = FALSE
  )
  rownames(name_map) <- original_names

  names(out)[names(out) == id] <- "id"
  names(out)[names(out) == treatment] <- "A"
  names(out)[names(out) == outcome] <- "Y"
  names(out)[names(out) == t_pros] <- "t_pros"
  names(out)[names(out) == os_time] <- "OS_time"
  names(out)[names(out) == death] <- "death"
  names(out)[names(out) == dropout] <- "dropout"

  for (i in seq_along(covariates)) {
    names(out)[names(out) == covariates[i]] <- paste0("X", i)
  }

  out$id <- as.character(out$id)
  out$A <- as.numeric(out$A)
  out$Y <- as.numeric(out$Y)
  out$t_pros <- as.numeric(out$t_pros)
  out$OS_time <- as.numeric(out$OS_time)
  out$death <- as.integer(out$death)
  out$dropout <- as.integer(out$dropout)

  x_cols <- grep("^X[0-9]+$", names(out), value = TRUE)
  for (x in x_cols) {
    out[[x]] <- as.numeric(out[[x]])
  }

  out <- out[!is.na(out$id) &
               !is.na(out$A) &
               !is.na(out$OS_time) &
               !is.na(out$death) &
               !is.na(out$dropout), ]

  if (any(out$OS_time < 0, na.rm = TRUE)) {
    stop("OS_time must be nonnegative.")
  }

  rownames(out) <- NULL

  return(list(data = out, name_map = name_map))
}
