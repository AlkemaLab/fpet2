#' Get smoothing standata
#'
#' @param fix_smoothing logical, whether to fix smoothing parameters
#' @param global_fit fit object from a global fit, used to extract smoothing parameters if fix_smoothing is TRUE
#' @param prefix prefix to add to the parameter names, defaults to ""
#'
#' @returns list with smoothing standata
#' @keywords internal
get_smoothing_standata <- function(fix_smoothing,
                                   global_fit,
                               prefix = ""){
  smoothing_data <- list()
  smoothing_data[[paste0(prefix, "Rho_fixed")]] <- numeric(0)
  smoothing_data[[paste0(prefix, "Tau_fixed")]] <- numeric(0)

  if (fix_smoothing) {
     if (is.null(global_fit)) {
       stop("fix_smoothing was set to TRUE, but a global_fit was not provided.")
      }
     # if (!smoothing) {
     #    stop("fix_smoothing was set to TRUE, but smoothing is FALSE.")
     # }
    # get estimates of the smoothing parameters rho and tau from the global fit
    smoothing_data[[paste0(prefix, "Rho_fixed")]] <-
      global_fit$post_summ %>%
      dplyr::filter(variable == paste0(prefix, "Rho[1]")) %>%
      pull(postmean)
    smoothing_data[[paste0(prefix, "Tau_fixed")]] <-
      global_fit$post_summ %>%
      dplyr::filter(variable == paste0(prefix, "Tau[1]")) %>%
      pull(postmean)
    # get estimates of the smoothing parameters rho and tau from the global fit
  }
  return(smoothing_data)
}

