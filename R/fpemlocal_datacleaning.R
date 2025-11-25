#re-using fpemlocal functions, copied from core data
#https://github.com/AlkemaLab/fpemlocal/blob/master/R/core_data.R
#note that impute_indicator() was related to imputing SEs and data_series_type_relabel() was related to getting indices for data series types.
# all internal function for now

fpemlocal_datacleaning <- function(dat){
  dat <-
    dat %>%
    ad_hoc_calculate_cp_trad() %>%
    # data_series_type_relabel() %>%
    # rounding up from zero etc
    dplyr::mutate(indicate_rounding_trad = indicate_rounding(contraceptive_use_traditional)) %>%
    dplyr::mutate(indicate_rounding_mod = indicate_rounding(contraceptive_use_modern)) %>%
    dplyr::mutate(contraceptive_use_traditional = round_from_zero(contraceptive_use_traditional)) %>%
    dplyr::mutate(contraceptive_use_modern = round_from_zero(contraceptive_use_modern)) %>%
    dplyr::mutate(contraceptive_use_any = round_from_zero(contraceptive_use_any)) %>%
    dplyr::mutate(unmet_need_any = round_from_zero(unmet_need_any)) %>%
    # recalculation after rounding
    ad_hoc_recalculate_cp_any() %>%
    ad_hoc_blankmodern_ifequals() %>%
    # creating a single column for subpopulation indicators
    dplyr::mutate(subpopulation_labels = subpopulation_labels(.)) %>%
    dplyr::mutate(subpopulation_descriptions = subpopulation_descriptions(.))
  return(dat)
  #levels(dat[["subpopulation_labels"]])  <- c("+", "-", "A", "F", "S-", "S+")
}




#' Round from zero
#'
#' Taken from UNPD/fpemlocal:
#' https://github.com/AlkemaLab/fpemlocal/blob/master/R/rounding_and_recalculation.R
#' @param x proportion to round
#' @return proportions less than 0.001 are rounded upwards to 0.001
#' @keywords internal
round_from_zero <- function(x) {
  min_round <- 0.001
  x[x < min_round] <- min_round
  return(x)
}


#' Indicate rounding
#'
#' Taken from UNPD/fpemlocal
#'
#' @param x proportion to check
#'
#' @returns indicator variable (1 if x < 0.001, else 0)
#' @export
#'
#' @keywords internal
indicate_rounding <- function(x) {
  min_round <- 0.001
  ind <- as.numeric(x < min_round)
  return(ind)
}


#' Ad hoc calculation of CP Traditional
#'
#' Taken from UNPD/fpemlocal: try to calculate cp traditional from any and modern
#'
#'
#' @param obs data frame with contraceptive_use_any, contraceptive_use_modern, contraceptive_use_traditional
#'
#' @returns data frame with contraceptive_use_traditional filled in where possible
#'
#' @keywords internal
ad_hoc_calculate_cp_trad <- function(obs) {
  obs <- obs %>%
    dplyr::mutate(contraceptive_use_traditional = ifelse(
      !is.na(contraceptive_use_any) & !is.na(contraceptive_use_modern) & is.na(contraceptive_use_traditional) & (contraceptive_use_modern != contraceptive_use_any),
      contraceptive_use_any - contraceptive_use_modern,
      contraceptive_use_traditional)
    )
  return(obs)
}

# re-calculated CP Any if CP Mod or CP Trad had been round up to make sure CP Any still equaled the sum of modern and traditional
ad_hoc_recalculate_cp_any <- function(obs) {
  obs <- obs %>%
    dplyr::mutate(recalc_cpany_indicator = (indicate_rounding_mod | indicate_rounding_trad) & !is.na(contraceptive_use_modern) & !is.na(contraceptive_use_traditional)) %>%
    dplyr::mutate(contraceptive_use_any = ifelse(recalc_cpany_indicator,
                                                 contraceptive_use_modern + contraceptive_use_traditional,
                                                 contraceptive_use_any
    )
    )
  return(obs)
}

#' Blank out modern use if modern cp and any CP are equal
#'
#' @param obs data frame with contraceptive_use_any, contraceptive_use_modern, contraceptive_use_traditional
#'
#' @returns data frame with contraceptive_use_modern set to NA if it equals contraceptive_use_any and contraceptive_use_traditional is NA
#'
#' @keywords internal
ad_hoc_blankmodern_ifequals <- function(obs) {
  obs <- obs %>%
    dplyr::mutate(
      contraceptive_use_modern = ifelse(
        contraceptive_use_any == contraceptive_use_modern &
          is.na(contraceptive_use_traditional) &
          !is.na(contraceptive_use_any) &
          !is.na(contraceptive_use_modern),
        NA,
        contraceptive_use_modern
      )
    )
  return(obs)
}



#' Assign labels for subpopulations
#'
#' Taken from UNPD/fpemlocal: see https://github.com/AlkemaLab/fpemlocal/blob/master/R/subpopulation_labels.R

#' @param contraceptive_use
#'
#' @returns labels for subpopulations
#' @keywords internal
subpopulation_labels <- function(contraceptive_use) {
  ifelse(
    contraceptive_use$age_group_bias == "+",
    "+",
    ifelse(
      contraceptive_use$age_group_bias == "-",
      "-",
      ifelse(
        contraceptive_use$age_group_bias == "?",
        "A",
        ifelse(
          contraceptive_use$has_traditional_method_bias == "Y",
          "F",
          ifelse(
            contraceptive_use$modern_method_bias == "-",
            "S-",
            ifelse(contraceptive_use$modern_method_bias == "+",
                   "S+",
                   "")
          )
        )
      )
    )
  )
}



#' Assign descriptions for subpopulations
#'
#' @param contraceptive_use data set with information on subgroups
#'
#' @returns descriptions for subpopulations
#'
#' @keywords internal
subpopulation_descriptions <- function(contraceptive_use) {
  ifelse(
    contraceptive_use$age_group_bias == "+",
    "+: age group bias",
    ifelse(
      contraceptive_use$age_group_bias == "-",
      "-: age group bias",
      ifelse(
        contraceptive_use$age_group_bias == "?",
        "A: age group bias",
        ifelse(
          contraceptive_use$has_traditional_method_bias == "Y",
          "F: traditional method bias",
          ifelse(
            contraceptive_use$modern_method_bias == "-",
            "S-: modern method bias",
            ifelse(contraceptive_use$modern_method_bias == "+",
                   "S+: modern method bias",
                   "")
          )
        )
      )
    )
  )
}
