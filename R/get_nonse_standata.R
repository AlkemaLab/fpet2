
#' Get stan data for non-sampling errors (nonse)
#'
#' @param fix_nonse logical, whether to fix nonse parameters to values from a global fit
#' @param global_fit fit object from a global fit, if fix_nonse is TRUE
#' @param source_index source_index data frame for the local fit
#' @param prefix prefix for nonse parameters, defaults to "DM1_"
#'
#' @returns list of stan data elements for nonse
#' @keywords internal
get_nonse_standata <- function(fix_nonse, global_fit,
                               source_index,
                               prefix = "DM1_"){
  nonse_data <- list()
  nonse_data[[paste0(prefix, "nonse_fixed")]] <-  numeric(0)
  parnames_outlier_hyper <- c("global_shrinkage_dm", "caux_dm", "sdbias",
                              "rho_pma")
  parnames_outlier_hyper_fixed <- paste0(parnames_outlier_hyper, "_fixed")
  # convoluted way to get each par in the list, with numeric(0) assigned
  outlier_fixed <- rep(0, length(parnames_outlier_hyper))
  names(outlier_fixed) <- parnames_outlier_hyper_fixed
  #nonse_data <- c(nonse_data, lapply(split(outlier_fixed,names(outlier_fixed)),unname))
  for (parname in parnames_outlier_hyper_fixed){
    nonse_data[[paste0(prefix, parname)]] <- numeric(0)
  }
  if (fix_nonse) {
    if (is.null(global_fit)) {
      stop("fix_nonse was set to TRUE, but a global_fit was not provided.")
    }
    if (!all(source_index$data_series_type %in% global_fit$source_index$data_series_type)) {
      stop("fix_nonse was set to TRUE, but data for local fit includes source types that were not in the data for the global fit.")
    }
    # get the nonse estimates from the global fit, ensuring that we get the
    # right ones, in the right order for sources in local fit
    for (parname in paste0(prefix, "nonse")){
      # to do: check ordering ok here?
      global_nonse_estimates <- global_fit$post_summ %>%
        dplyr::filter(variable_no_index == parname) %>%
        mutate(data_series_type = global_fit$source_index$data_series_type)
      nonse_data[[paste0(parname, "_fixed")]] <-
        source_index |>
        dplyr::left_join(global_nonse_estimates, by = "data_series_type") |>
        dplyr::pull(postmean)
    }
    # for outliers
    parnames_outlier_hyper_tofix <-
      expand.grid(
        prefix_use = prefix,
        param = parnames_outlier_hyper) %>%
        dplyr::mutate(prefixed_param_name = paste0(prefix_use, param)) %>%
        pull(prefixed_param_name)
    for (parname in parnames_outlier_hyper_tofix){
      nonse_data[[paste0(parname, "_fixed")]] <-
        global_fit$post_summ %>%
        dplyr::filter(variable == paste0(parname, "[1]")) %>%
        pull(postmean)
    }
  }# end fixing dm pars
  return(nonse_data)
}
