#' Simulated terminal trend data for `TTM()`
#'
#' A simulated long-format data set with repeated measurements of a
#' longitudinal outcome before death, dropout, or censoring. It is used in the
#' examples and the vignette of [TTM()] and [TTM_linear()]. The script used to
#' generate the data is in the `data-raw` folder of the source repository.
#'
#' @format A data frame with 3681 rows (1000 subjects) and 17 variables:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{A}{Binary treatment indicator (0 = control, 1 = treatment).}
#'   \item{X1}{Binary baseline covariate.}
#'   \item{X2}{Non-negative continuous baseline covariate.}
#'   \item{u}{Subject-level random effect generated in the simulation.}
#'   \item{b_longi}{Random intercept for the longitudinal outcome.}
#'   \item{randomslope}{Random slope generated for the longitudinal outcome.}
#'   \item{OS_time}{Observed time, the minimum of the death, dropout, and
#'     censoring times (months).}
#'   \item{event}{Event indicator: 0 = censored, 1 = death, 2 = dropout.}
#'   \item{death}{Death indicator (1 = death observed).}
#'   \item{dropout}{Dropout indicator (1 = dropout observed).}
#'   \item{death_time}{Simulated (latent) death time.}
#'   \item{dropout_time}{Simulated (latent) dropout time.}
#'   \item{censoring_time}{Simulated (latent) administrative censoring time.}
#'   \item{t_pros}{Prospective measurement time since baseline.}
#'   \item{t_retro}{True retrospective time, `death_time - t_pros`.}
#'   \item{Y}{Longitudinal outcome.}
#' }
#' @source Simulated data.
#' @examples
#' data("TTM_data", package = "TTM")
#' head(TTM_data)
#' table(TTM_data$event[!duplicated(TTM_data$id)])
"TTM_data"

#' Simulated retrospective joint modeling data for `RetroJM()`
#'
#' A simulated long-format data set with repeated measurements of a
#' longitudinal outcome before death or censoring. It is used in the
#' vignette of [RetroJM()]. It has the same structure as [TTM_data].
#'
#' @format A data frame with 3586 rows (1000 subjects) and 17 variables:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{A}{Binary treatment indicator (0 = control, 1 = treatment).}
#'   \item{X1}{Binary baseline covariate.}
#'   \item{X2}{Non-negative continuous baseline covariate.}
#'   \item{u}{Subject-level random effect generated in the simulation.}
#'   \item{b_longi}{Random intercept for the longitudinal outcome.}
#'   \item{randomslope}{Random slope generated for the longitudinal outcome.}
#'   \item{OS_time}{Observed time, the minimum of the death, dropout, and
#'     censoring times (months).}
#'   \item{event}{Event indicator: 0 = censored, 1 = death, 2 = dropout.}
#'   \item{death}{Death indicator (1 = death observed).}
#'   \item{dropout}{Dropout indicator (1 = dropout observed).}
#'   \item{death_time}{Simulated (latent) death time.}
#'   \item{dropout_time}{Simulated (latent) dropout time.}
#'   \item{censoring_time}{Simulated (latent) administrative censoring time.}
#'   \item{t_pros}{Prospective measurement time since baseline.}
#'   \item{t_retro}{True retrospective time, `death_time - t_pros`.}
#'   \item{Y}{Longitudinal outcome.}
#' }
#' @source Simulated data.
#' @examples
#' data("Retro_data", package = "TTM")
#' head(Retro_data)
"Retro_data"
