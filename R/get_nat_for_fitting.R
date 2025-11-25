
#' Get national data for subnational fitting
#'
#' @param national_dat_df National survey data
#' @param survey_df_use Subnational survey data used in fitting
#'
#' @returns Dataframe with national data not in subnational survey data
#' @keywords internal
get_nat_for_fitting <- function(national_dat_df, survey_df_use) {
  # call by marital group
  # returns national data that is not in the survey data
  dat_tofit <- national_dat_df %>%
    # filter out PMA from national data to consider
    dplyr::filter(data_series_type != c("PMA")) %>%
    # filter out start_end year and data_source type that are already in survey by region
    dplyr::anti_join(survey_df_use, by = c("division_numeric_code", "start_date", "end_date", "data_series_type"))
  return(dat_tofit)
}
