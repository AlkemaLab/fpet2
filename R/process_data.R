
#' Process survey data
#'
#' @param dat Survey data
#' @param regions_dat File with region meta data
#' @param start_year Filter out data prior to this year
#' @param is_married Logical, select married women data?
#'
#' @returns Processed survey data
#' @export
#'
process_data <- function(dat,
                         regions_dat,
                         start_year = 1970,
                         is_married = TRUE

  ){
  # data processing to prepare data for model fitting

  ### filtering for marital status and start year
  dat <- dat %>%
    dplyr::filter(is_in_union == ifelse(is_married, "Y", "N")) %>%
    dplyr::mutate(year = floor((start_date + end_date) / 2)) %>%
    dplyr::filter(year >= start_year)

  if (dim(dat)[1] == 0) {
    print("No data")
    return(NULL)
  }

  # age-range here is column that indicates whether data to be included in this model
  # (actual age ranges are used to define biases)
  dat <- dat %>%
    dplyr::filter(age_range == "15-49")

  # set unmet need data less than .1% as NA
  dat <- dat %>%
    mutate(unmet_need_any = ifelse(unmet_need_any < 0.001, NA, unmet_need_any))

  # apply UNPD processing
  dat <- fpemlocal_datacleaning(dat)

  # calculate things
  dat <- dat %>%
    mutate(demand = contraceptive_use_any + unmet_need_any,
           # default is unmet=. unmet modern, and unmet_any
           unmet = ifelse(!is.na(unmet_need_modern), unmet_need_modern, unmet_need_any + contraceptive_use_traditional),
           unmet_any = unmet - contraceptive_use_traditional,
           use_any = contraceptive_use_modern + contraceptive_use_traditional) %>%
    mutate(unmet_overnotmodern = unmet/(1-contraceptive_use_modern),
           trad_overnotmodern = contraceptive_use_traditional/(1-contraceptive_use_modern),
           trad_overunmet = contraceptive_use_traditional/unmet,
           unmetany_overnonuse = unmet_any/(1-use_any),
           demand_satisfied_modern = contraceptive_use_modern/demand,
           demand_satisfied_any = contraceptive_use_any/demand,
           ratio_modern_any = contraceptive_use_modern/contraceptive_use_any)

  # adding in country names and cluster (region) info
  dat <-
    dat %>%
    left_join(regions_dat, by = "division_numeric_code")
  if (any(is.na(dat$cluster))){
    print("Some clusters are NA, these observations are excluded for now.")
    print("Observations with missing clusters")
    print(dat %>% dplyr::filter(is.na(cluster)) %>% dplyr::select(iso, year, cluster))
    dat <- dat %>% dplyr::filter(!is.na(cluster))
  }
  if (!("record_id_fixed" %in% names(dat))){
    dat$record_id_fixed <- paste0("TMPID_", seq(1, dim(dat)[1]))
  }


  dat <- get_se_sequentialdm(dat)
  dat <- add_bias_info(dat, married = is_married)
  dat <- assign_outliers(dat, married = is_married)
  # for running the model, we use this function to assign dummies
  # (can be moved inside run_step)
  # dat2 <- assign_outliers(dat2)
  #  dat2 %>% dplyr::select(iso, year, data_series_type, possible_outlier, possible_outlier_userinput, nooutlier, any_bias, record_id_fixed)




  # TO DO: come back to check PMA data, needed again in fit_model?
  # check that pma doesn't have missing trad and unmet_modern entries for ALL steps
  tmp <- dat %>%
    dplyr::filter(data_series_type == "PMA") %>%
    dplyr::filter(!is.na(contraceptive_use_traditional) & is.na(unmet))
  if (dim(tmp)[1] != 0){
    stop("PMA data on trad is given but PMA data on unmet is missing, we don't allow for that.")
  }



  return(dat%>%
           mutate(indic = "fp"))
}

