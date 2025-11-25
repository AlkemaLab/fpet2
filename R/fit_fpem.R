
#' Fit FPET national or subnational local models for ALL women, with or without EMUs.
#'
#' @param survey_df Survey data for model fitting, filtered to pop of interest
#' @param population_df Population data, filtered to pop of interest
#' @param service_statistic_df EMU data, filtered to pop of interest
#' @param subnational logical, TRUE if fitting to subnational data.
#' For subnational fits, survey_df may contain all or just 1 region.
#' @param national_dat_df National data to be used/displayed for subnational runs,
#' filtered to pop of interest. Model is fitted to surveys in national_dat_df that are not in survey_df.
#' @param chains Number of chains for Stan
#' @param iter_sampling Number of iterations for sampling
#' @param iter_warmup Number of iterations for warmup
#' @param seed Random seed for Stan
#' @param regions_dat Data frame with region information
#' @param end_year End year for the model, default is 2030
#'
#' @returns results list with estimates, samples, data_allwomen, observations,
#'  and dat_emu = dat_emu
#' @export
#'
fit_fpem  <- function(
    survey_df,
    population_df,
    service_statistic_df = NULL, # filtered to division numeric code already
    subnational = FALSE,
    national_dat_df = NULL,
    chains = 4,
    iter_sampling = 200,
    iter_warmup = 150,
    seed = 1234,
    regions_dat,
    end_year = 2030
) {

  print("This is the FPET2025 version 1.2 (released November 24, 2025)!")
  if (subnational) {
    # we need to get geo_units in same order
    survey_df <- survey_df %>%
      dplyr::arrange(region_code)
  }
  if (sum(is.na(survey_df$contraceptive_use_modern)) > 0){
    print("We remove observations that do not have mCPR observed.")
    survey_df <- survey_df %>%
      dplyr::filter(!is.na(contraceptive_use_modern))
  }

  # create list with the data inputs for married and unmarried
  data_allwomen <- list()
  nodata_unmarried <- FALSE
  for (is_married in c(TRUE, FALSE)){
    marital_status <- ifelse(is_married, "married", "unmarried")
    survey_df_use <- process_data(survey_df, regions_dat,
                                  is_married = is_married)
    if (is.null(survey_df_use) & marital_status == "unmarried") {
      nodata_unmarried <- TRUE
      print("No data for unmarried women, we use married women data and impute NAs")
      survey_df_use <- survey_df %>%
                 mutate(is_in_union  = "N",
                        # need to make sure that the data type was included
                        data_series_type = "DHS",
               contraceptive_use_modern = NA,
               contraceptive_use_traditional = NA,
               contraceptive_use_any = NA,
               unmet_need_modern = NA,
               unmet_need_any = NA)
      survey_df_use <- process_data(survey_df_use, regions_dat,
                                    is_married = is_married)
    }
    # get population_df_use per marital group, needed only for subnational
    if (subnational) {
      population_df_use <- population_df %>%
        dplyr::filter(is_in_union == ifelse(is_married, "Y", "N"))
    } else {
      population_df_use <- NULL
    }

    runstep <- ifelse(subnational, "local_subnational", "local_national")
    # fit <- readRDS(file = here::here("data-raw/internal",
    #                                paste0(
    #                                  ifelse(subnational, "fit3_fromnat_", "fit2_"),
    #                                  marital_status, ".rds")))
    fitname <- paste0(ifelse(subnational, "fit3_fromnat_", "fit2_"),
                      marital_status)
    assign("fit", get(fitname))
    if (subnational){
      survey_df_use <- survey_df_use %>%
        mutate(is_aggr = FALSE)
      if (!is.null(national_dat_df)){
        if (sum(is.na(national_dat_df$contraceptive_use_modern)) > 0){
          print("We remove national observations that do not have mCPR observed.")
          national_dat_df <- national_dat_df %>%
            dplyr::filter(!is.na(contraceptive_use_modern))
        }
        national_dat_df_use <- process_data(national_dat_df, regions_dat,
                                      is_married = is_married)
        if (!is.null(national_dat_df_use)){
          national_dat_forfitting <- get_nat_for_fitting(national_dat_df_use, survey_df_use)
          if (dim(national_dat_forfitting)[1] > 0){
            # move from nat data to survey data and add isagrr column
            national_dat_df <- national_dat_df %>%
              anti_join(national_dat_forfitting, by = c("start_date",
                                                         "end_date",
                                                         "is_in_union"))
            survey_df_use <- bind_rows(
              survey_df_use,
              national_dat_forfitting %>% mutate(is_aggr = TRUE,
                                                 # region_code is set to NA in fit_model such that not included in geo_index
                                                 region_code = "National")
            )
          }
        }
      }
    }   # end adding aggregate column and considering national data for fitting

    data_allwomen[[marital_status]] <-  fit_model(
      iso_select  = survey_df_use$iso[1],
      runstep = runstep,
      area = ifelse(subnational, "region_code", "iso"),
      # tmp testing
      #add_aggregates = FALSE,
      add_aggregates = TRUE, # not used for national
      population_df = population_df_use,
      global_fit = fit,
      survey_df = survey_df_use,
      is_married = is_married,
      add_sample = FALSE,
      end_year = end_year
    )
  } # end processing by marital group

  if (is.null(service_statistic_df)){
    dat_emu <- NULL
    emu_list <- NULL

  } else {
    #hyper_param <- readRDS("data-raw/internal/emu_global_hyperparameters.RDS")
    #assign("hyper_param", emu_global_hyperparameters)
    hyper_param <- fpet2::emu_global_hyperparameters
    combined_list <- service_statistics_getstandata(
      service_statistic_df =
        service_statistic_df %>% left_join(
          regions_dat %>% dplyr::select(division_numeric_code, iso)),
      time_index = data_allwomen[["married"]]$time_index, # to get time index
      hyper_param = hyper_param,
      geo_unit_index = data_allwomen[["married"]]$geo_unit_index %>%
        dplyr::select(any_of(c("iso", "region_code", "c")))
    )
    dat_emu <- combined_list$dat_emu
    emu_list <- combined_list$emu_list
  }

  fit_all <- fit_model_allwomen(
                      runstep = runstep,
                      data_allwomen = data_allwomen,
                      population_df = population_df,
                      emu_list = emu_list,
                      end_year = end_year,
                      seed = seed
                    )
  #return(fit_all)
  print("Chains finished, now calculating estimates (can take a little while)")
  results <- get_estimates_fpem(fit_all$samples,
                            geo_unit_index = data_allwomen[[1]]$geo_unit_index,
                            time_index = data_allwomen[[1]]$time_index,
                            subnational = subnational,
                            return_samples = TRUE)#FALSE)
  if (!nodata_unmarried){
    observation_df <- bind_rows(
      add_uncertainty_in_obs_allwomen(fit_all$samples,
                                      data_allwomen,
                                      is_married = TRUE),
      add_uncertainty_in_obs_allwomen(fit_all$samples,
                                      data_allwomen,
                                      is_married = FALSE))
  } else {
    observation_df <- add_uncertainty_in_obs_allwomen(fit_all$samples,
                                      data_allwomen,
                                      is_married = TRUE)
  }
  if (subnational){
    # if there was aggregate data included in the fitting, region_code is NA
    observation_df$region_code <- ifelse(is.na(observation_df$region_code),
                                          "National",
                                          observation_df$region_code)
    # if there is additional data not used in the fitting, add for plotting
    if (!is.null(national_dat_df)){
      print("Processing of additional national data not used in fitting, to add to observation_df and show in plots.")
      if (dim(national_dat_df %>% dplyr::filter(is_in_union == "Y"))[1] > 0){
        observation_df  <- observation_df  %>%
          bind_rows(add_natdata_for_aggregates(national_dat_df, regions_dat, is_married = TRUE))
      }
      if (dim(national_dat_df %>% dplyr::filter(is_in_union == "N"))[1] > 0){
        observation_df  <- observation_df  %>%
          bind_rows(add_natdata_for_aggregates(national_dat_df, regions_dat, is_married = FALSE))
      }
    }
  }
  # add population count

  # first add all women to the pop count object
  group_col <- c("division_numeric_code", "year")
  if (subnational){
    group_col <- c(group_col, "region_code")
  }
  population_df <-
    bind_rows(
      population_df,
      population_df |>
        dplyr::filter(is_in_union %in% c("N", "Y")) |>
        group_by(across(all_of(group_col))) |>
        summarize(
          is_in_union = "all",
          population_count = sum(population_count)
        )
    )

  if (subnational){
    # also add pop numbers for national
    population_df <-
      bind_rows(
        population_df,
        population_df %>%
          group_by(division_numeric_code, year, is_in_union) %>%
          summarize(population_count = sum(population_count)) %>%
          mutate(region_code = "National")
      )
  }

  results$estimates <-
    results$estimates %>%
    rename(proportion = value) %>%
    dplyr::select(-measure) %>%
    left_join(population_df %>%
                mutate(marital_status = case_when(
                  is_in_union == "all" ~ "all",
                  is_in_union == "Y" ~ "married",
                  is_in_union == "N" ~ "unmarried"
                ))) %>%
    # add counts for all except ratios
    mutate(population_count =
             case_when(
               indicator %in% c("demand_satisfied" , "demand_satisfied_modern" , "ratio_modern_any") ~ NA,
               TRUE ~ proportion*population_count)) %>%
    dplyr::select( - is_in_union) %>%
    tidyr::pivot_longer(cols = c(proportion, population_count), names_to = "measure") %>%
    dplyr::mutate(division_numeric_code = survey_df$division_numeric_code[1])
  results <- c(
    results,
    list(
      samples = fit_all$samples,
      data_allwomen = data_allwomen,
     observations = observation_df,
    # geo_unit_index = data_allwomen[[1]]$geo_unit_index, # if we need it?
     dat_emu = dat_emu
   ))
  return(results)
}


