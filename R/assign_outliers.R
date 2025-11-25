


#' Assign outliers
#'
#' Function to assign outliers based on a set of rules.
#' If columns possible_outlier or possible_outlier_userinput are present, these are used to define outliers.
#'
#' @param fp_dat survey data
#' @param outlier_record_ids record_ids of outliers (if any)
#' @param married boolean, needed when bias needs to be assessed
#'
#' @returns data frame with new column nooutlier (1 = not an outlier, 0 = outlier)
#' @keywords internal
#'
assign_outliers <- function(fp_dat,
                            outlier_record_ids = NULL, # needed for 1b, when adding outlier candidates
                            married # needed when bias needs to be assessed
                            ){
  nat_data <- fp_dat

  if ("possible_outlier" %in% names(nat_data)){
    if ("possible_outlier_userinput" %in% names(nat_data)){
      if (any(!is.na(nat_data$possible_outlier_userinput))){
        print("We define possible outliers based on columns possible_outlier and possible_outlier_userinput")
        nat_data <- nat_data %>%
          mutate(nooutlier = 1 - ifelse(is.na(possible_outlier_userinput), possible_outlier, possible_outlier_userinput))
        return(nat_data)
      }
    }
    if (any(is.na(nat_data$possible_outlier))){
        print("Column possible outliers has NAs, it will be regenerated")
    } else {
      print("We define possible outliers based on column possible_outlier")
      nat_data <- nat_data %>%
          mutate(nooutlier = 1 - possible_outlier)
      return(nat_data)
    }
  } else {
    print("Column possible_outlier not included in data, no pre-defined outliers used")
  }


  # from earlier global tests
  if (!("outlier_la" %in% names(nat_data))){
    #print("Column outlier_la not included in data, no pre-defined outliers used")
    nat_data <- nat_data %>%
      mutate(outlier_la = NA)
  }

  # need any_bias for classification
  if (!("any_bias" %in% names(nat_data))){
    nat_data <- add_bias_info(nat_data, married = married)
  }

  nat_data <- nat_data %>%
    mutate(isafter1990 = ifelse(start_date > 1990, 1, 0),
           yesoutlier = case_when(
             outlier_la ==1 ~ 1,
             any_bias == 1 ~ 1,
             data_series_type == "DHS" & !isafter1990 ~ 1,
             .default = 0
           )) %>%
#    group_by(name_country) %>%
    group_by(division_numeric_code) %>%
    mutate(ndhs_nonoutlying = sum(data_series_type == "DHS" & !yesoutlier),
           # apply recency by counting after 1990
           nnationalsurvey_nonoutlying = sum(data_series_type == "National survey" & !yesoutlier & isafter1990),
           nother_nonoutlying = sum(data_series_type == "Other"  & !yesoutlier & isafter1990)
    ) %>%
    ungroup() %>%
    mutate(
      nooutlier = case_when(
        yesoutlier == 1  ~ 0,
        data_series_type == "DHS"  ~ 1,
        # do assign as reference also prior to 1990
        ndhs_nonoutlying == 0 & nnationalsurvey_nonoutlying >= nother_nonoutlying &
          data_series_type == "National survey" ~ 1,
        ndhs_nonoutlying == 0 & nnationalsurvey_nonoutlying < nother_nonoutlying &
          data_series_type == "Other" ~ 1,
        .default = 0
      )) %>%
    group_by(division_numeric_code) %>%
#    group_by(name_country) %>%
    # if nothing else, use the most recent non-outlying point in a country
    mutate(noreferenceyet = ifelse(sum(nooutlier) == 0 , 1, 0),
           most_recent_notoutlying = max(c(0, # just a small number added to avoid -Inf warnings
                                            start_date[yesoutlier == 0]))) %>%
    ungroup() %>%
    mutate(nooutlier = ifelse(noreferenceyet & start_date == most_recent_notoutlying, 1, nooutlier))

  if (!is.null(outlier_record_ids)){
    print("we use old outlier info and outlier_record_ids")
    if (!("record_id_fixed" %in% names(nat_data))){
      stop("Column record_id_fixed not included in data, we need that!")
    }
    if (any(is.na(nat_data$record_id_fixed))){
      stop("Missing values in record_id_fixed not allowed!")
    }
    #print(mean(nat_data$nooutlier))
    nat_data <- nat_data %>%
      mutate(nooutlier = ifelse (record_id_fixed %in% outlier_record_ids | !nooutlier, 0, 1))
  }
  nat_data
}


#' Add possible-outlier related columns
#'
#' To survey data, add columns possible_outliers and possible_outlier_userinput. Possible outliers are assessed
#' using the `assign_outliers` function.
#'
#' @param dat survey data
#'
#' @returns data frame with new columns possible_outlier (1 = possible outlier, 0 = not a possible outlier) and possible_outlier_userinput (NA)
#' @keywords internal
add_outlier_related_columns <- function(dat){
  dplyr::bind_rows(
    dat %>%
      dplyr::filter(is_in_union == "Y") %>%
      mutate(possible_outlier = 1-  assign_outliers( ., married = TRUE) %>% pull(nooutlier)),
    dat %>%
      dplyr::filter(is_in_union == "N") %>%
      mutate(possible_outlier = 1 - assign_outliers( ., married = FALSE) %>% pull(nooutlier))
  ) %>%
    mutate(possible_outlier_userinput = NA)
}


#' Use global fit 1a residuals to update assignment of outlier errors
#'
#' In this function, we calculate residuals from global fit 1a to
#' update the `nooutlier` column for use in further model fitting.
#' In this function, `nooutlier` is set to 0 for observations that are outlying in 1a,
#' here defined as observations with absolute standardized residuals above the 90th percentile.
#' By setting `nooutlier` to zero, outlier errors are assigned.
#' (Naming of the column can use an update!)
#'
#'
#' @param fit global fit 1a
#' @param is_married logical
#'
#' @returns data frame with updated nooutlier column.
#' @export
#'
update_nooutlier_in_global_fit <- function(fit, is_married){
  # takes in a fit, outputs data from that fit with outliers assigned

  resid <- get_residuals(fit, is_married)
  cutoff_y <- quantile(abs(resid$mean_st_res_modern), 0.9, na.rm = TRUE)
  cutoff_unmet <- quantile(abs(resid$mean_st_res_unmetovernonmodern), 0.9, na.rm = TRUE)
  resid2 <- resid %>%
    mutate(isoutlier = ifelse(abs(mean_st_res_modern) > cutoff_y |
                                abs(mean_st_res_unmetovernonmodern) > cutoff_unmet,
                              1, 0)) %>%
    mutate(isoutlier = ifelse(!is.na(isoutlier), isoutlier, 0))

  dat <- fit$original_data


  # dat$nooutlier_updated <- ifelse(resid2$isoutlier == 1, 0, dat$nooutlier)
  # dat %>%
  #   dplyr::filter(nooutlier != nooutlier_updated) %>%
  #   dplyr::select(iso, data_series_type, year, nooutlier, nooutlier_updated)
  dat <- dat %>%
    dplyr::mutate(nooutlier = ifelse(resid2$isoutlier == 1, 0, nooutlier)) %>%
    # in 1a, we set dhs to dhs0 so revert that here if needed
    dplyr::mutate(data_series_type = ifelse(data_series_type == "DHS0", "DHS", data_series_type))
  dat <- dat %>%
    dplyr::select(all_of(names(fit$original_data)))
  return(dat)
}


