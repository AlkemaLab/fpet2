
#' Add national data for aggregates
#'
#' This function processes national data and adds it to the aggregate dataset.
#'
#' @param national_dat_df national survey data prior to processing
#' @param regions_dat region metadata
#' @param is_married boolean
#'
#' @returns Processed national data with region_code set to "National"
#' @keywords internal
add_natdata_for_aggregates <- function(national_dat_df,
                                       regions_dat,
                                       is_married){
  nat_data_use <- process_data(national_dat_df,
                               regions_dat,
                               is_married = is_married) %>%
    mutate(region_code = "National") %>%
    dplyr::mutate(
      est_contraceptive_use_modern = contraceptive_use_modern,
      est_contraceptive_use_traditional = contraceptive_use_traditional,
      est_contraceptive_use_any = contraceptive_use_any,
      est_unmet_need_any = unmet_need_any,
      est_demand = demand,
      est_demand_satisfied_modern = demand_satisfied_modern,
      est_demand_satisfied_any = demand_satisfied_any,
      est_ratio_modern_any = ratio_modern_any,
      est_unmet_need_modern = unmet # unmet is unmet_modern in data
    )
  return(nat_data_use)
}
