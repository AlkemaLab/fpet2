
#' Get Stan data for service statistics data
#'
#' This function is called from fit_fpem, to prepare service statistics data for use in model fitting
#'
#' @param service_statistic_df emus filtered to pop (eg country)
#' @param time_index time_index used in model fit
#' @param hyper_param hyper parameters for EMUs to be used in model fit
#' @param geo_unit_index geo_unit_index used in model fit
#'
#' @returns list with dat_emu (for plotting) and emu_list (to pass to stan)
#' @keywords internal
service_statistics_getstandata <- function(service_statistic_df, # filtered to pop (eg country)
                                           time_index, # to get time index
                                           hyper_param,
                                           geo_unit_index
){
  dat_emu <-
    service_statistic_df %>%
    dplyr::mutate(emu_lower = pmax(0, emu - qnorm(0.975)*sd_emu),
                  emu_upper = pmin(1, emu + qnorm(0.975)*sd_emu),
                  sd2_emu_roc = sd_emu_roc^2) %>%
    left_join(time_index, by = "year") %>%
    mutate(emu_for_allwomen_j = ifelse(pop_type %in% c("MW","Married women"), FALSE, TRUE)) %>%
    # exclude data prior to the most recent survey? right now we don't
    rename(t_emu_j = t) %>%
    left_join(geo_unit_index) %>%
    rename(c_emu_j = c)
  if (dim(dat_emu)[1] == 0){
    dat_emu <- NULL
    emu_list <- NULL
  } else {
    if (anyNA(dat_emu$c_emu_j)){
      stop("some mismatch between regions provided in EMU data vs those in survey data.")
    }
    dat_emu_countrytype <-
      dat_emu %>%
      dplyr::select(ss_type, division_numeric_code) %>%
      dplyr::distinct() %>%
      dplyr::left_join(hyper_param) %>%
      dplyr::mutate(ss_type_index = row_number()) %>%
      dplyr::select(ss_type_index, ss_type, mean_log_sigma_type, hierarchical_sigma)

    dat_model <- dat_emu %>%
      dplyr::filter(!is.na(emu_roc))  %>% # we don't filter yet for dat_emu as we do want start years in the plot
      dplyr::left_join(dat_emu_countrytype) %>%
      dplyr::rename(ss_type_index_j = ss_type_index) %>%
      dplyr::select(emu_roc, sd2_emu_roc, emu_for_allwomen_j, c_emu_j, t_emu_j, ss_type_index_j)

    emu_list <- c(list(N_emu = dim(dat_model)[1]),
                  as.list(dat_model),
                  Ncountrytype = dim(dat_emu_countrytype)[1],
                  as.list(dat_emu_countrytype %>%
                            dplyr::select(mean_log_sigma_type, hierarchical_sigma)))
  }

  return(list(dat_emu = dat_emu, emu_list = emu_list))
}
