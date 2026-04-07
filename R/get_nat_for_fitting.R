
#' Get national data for subnational fitting
#'
#' @param national_dat_df National survey data, for specific marital group and div code
#' @param survey_df_use Subnational survey data used in fitting, for specific marital group and div code
#'
#' @returns Dataframe with national data if not present for all regions in subnational survey data
#' @keywords internal
get_nat_for_fitting <- function(national_dat_df, survey_df_use) {
  # first summarize the number of observations of a specific survey in survey data
  surveys_inallregionsinsubnat <-
    survey_df_use %>%
    # summarize the number of regions that have each specific survey
    # called by marital group and country, so unique follows from period and series type
    group_by(start_date, end_date, data_series_type) %>%
    summarise(n_region = n()) %>%
    ungroup() %>%
    # then keep just the surveys that are present in all regions
    filter(n_region == length(unique(survey_df_use$region_code)))

  dat_tofit <- national_dat_df %>%
    # filter out PMA from national data to consider
    dplyr::filter(data_series_type != c("PMA")) %>%
    ## filter out surveys that are already in survey data by region
    ## dplyr::anti_join(survey_df_use, by = c("division_numeric_code", "start_date", "end_date", "data_series_type"))
    # update 2026/4/7: only filter out national data if it's present in all regions
    dplyr::anti_join(surveys_inallregionsinsubnat,
                     by = c("start_date", "end_date", "data_series_type"))
  return(dat_tofit)
}
