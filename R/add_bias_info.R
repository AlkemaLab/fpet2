#' Add bias info
#'
#' Add columns for any_bias, trad_bias, tradorany_bias
#'
#' @param data survey data
#' @param married boolean
#'
#' @returns tibble based on data with added columns for any_bias, trad_bias, tradorany_bias
#'
#' @keywords internal
add_bias_info <- function(data, married){
  data %>%
    # added for 2024 data set
    mutate(modern_method_bias = ifelse(is.na(modern_method_bias), "None", modern_method_bias)) %>%
    mutate(any_bias = ifelse(modern_method_bias == "None" # modern applies to all
                             & has_geographical_region_bias == "N"
                             & has_non_pregnant_and_other_positive_biases == "N"
                             & age_group_bias == "None"
                             #(TRUE|FALSE)&TRUE
                             & ( (married == TRUE & group_type_relative_to_baseline == "MW") |
                                   (married  == FALSE & group_type_relative_to_baseline == "UW") ),
                             0, 1),
           trad_bias = ifelse(has_traditional_method_bias == "N"
                              & has_absence_of_probing_questions_bias == "N", 0, 1),
           tradorany_bias = ifelse(any_bias ==1 | trad_bias ==1, 1, 0))
}
