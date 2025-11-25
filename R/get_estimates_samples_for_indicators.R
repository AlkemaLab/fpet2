

#' Get samples for indicators
#'
#' This function extracts samples for specified indicators from a posterior samples object
#'
#' @param indicators vector of indicator names (e.g., c("contraceptive_use_modern", "unmet_need_modern"))
#' @param has_geo logical, if TRUE indicators are [geo_unit, time] indicators; if FALSE, just [time] (for aggregates)
#' @param samples posterior samples object
#' @param geo_unit_index data frame with geo unit information (required if has_geo is TRUE)
#' @param time_index data frame with time information
#' @param stringtoremove string to remove from indicator names
#' @param rename_modern logical, if TRUE renames indicators "modern", "unmet", "trad" to full names
#'
#' @returns tibble with samples for the specified indicators
#' @keywords internal
get_samples_for_indicators <- function(indicators,
                                       has_geo = TRUE, # if TRUE indicators are [geo_unit, time] indicators
                                       # if false, just [time] (for aggregates)
                                       samples,
                                       geo_unit_index = NULL,
                                       time_index,
                                       stringtoremove = "",
                                       rename_modern = FALSE
                                       ){
  if (has_geo){
    parameter_pattern <- paste0("^(",  paste(indicators, collapse = "|"),  ")\\[(\\d+),(\\d+)\\]$")
    samples_all <-
      samples$draws(indicators) |>
      posterior::as_draws_df() |>
      tidyr::pivot_longer(cols = !(.chain:.draw)) |>
      dplyr::rename_with(~ stringr::str_replace(.x, "^\\.", "")) |>
      dplyr::mutate(
        stringr::str_match(name, parameter_pattern) |>
          as_tibble(.name_repair = ~ c("matched", "parameter", "c", "t"))
      ) |>
      dplyr::mutate(
        c = as.integer(c),
        t = as.integer(t)
      ) |>
      dplyr::inner_join(geo_unit_index %>% select(- any_of(c("is_unmarried_sexual_activity", "hier_regional_unmarried",
                                                             "cluster", "cluster_unmarried", "subcluster"))), by = "c") |>
      dplyr::inner_join(time_index, by = "t")
  } else {
    parameter_pattern <- paste0("^(",  paste(indicators, collapse = "|"),  ")\\[(\\d+)\\]$")
    samples_all <-
      samples$draws(indicators) |>
      posterior::as_draws_df() |>
      tidyr::pivot_longer(cols = !(.chain:.draw)) |>
      dplyr::rename_with(~ stringr::str_replace(.x, "^\\.", "")) |>
      dplyr::mutate(
        stringr::str_match(name, parameter_pattern) |>
          as_tibble(.name_repair = ~ c("matched", "parameter",  "t"))
      ) |>
      dplyr::mutate(
        t = as.integer(t)
      ) |>
      dplyr::inner_join(time_index, by = "t")
  }
  samples_all <- samples_all %>%
    dplyr::select(-chain, -iteration, - name, -matched) %>%
    dplyr::rename(indicator = parameter) %>%
    # remove agg
    dplyr::mutate(indicator = gsub("_aggr", "", indicator)) %>%
    dplyr::mutate(indicator = gsub(stringtoremove, "", indicator))

  if (rename_modern){
    samples_all <- samples_all %>%
    tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
    tidyr::pivot_longer(cols = c(modern, unmet, trad), names_to = "indicator", values_to = "value") %>%
      dplyr::mutate(indicator = case_when(
      indicator == "modern" ~ "contraceptive_use_modern",
      indicator == "unmet" ~ "unmet_need_modern",
      indicator == "trad" ~ "contraceptive_use_traditional",
      TRUE ~ indicator))
  }
  return(samples_all)
}


#' Summarize family planning indicators helper
#'
#' This function summarizes family planning indicators from samples
#'
#' @param samples  tibble with columns indicator and value, where
#' indicator includes (unmet_need_any, contraceptive_use_modern,
#' contraceptive_use_traditional)
#'
#' @returns tibble with additional family planning indicators calculated
#' @keywords internal
summarize_fp_helper <- function(samples){
  ncolumns <- dim(samples)[2]
  samples_all <- samples %>%
    tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
    dplyr::select(-draw) %>%
    mutate(contraceptive_use_any = contraceptive_use_modern + contraceptive_use_traditional,
           non_use = 1 - contraceptive_use_any,
           unmet_need_any = unmet_need_modern - contraceptive_use_traditional,
            # unmet_need_modern = unmet_need_any + contraceptive_use_traditional,
           demand = contraceptive_use_any + unmet_need_any,
           # demand_modern added here for consistency with fpet1, note that it is just demand so could be removed
           demand_modern = demand,
           demand_satisfied = contraceptive_use_any/demand_modern,
           demand_satisfied_modern = contraceptive_use_modern/demand_modern,
           no_need = 1 - demand,
           ratio_modern_any = contraceptive_use_modern/contraceptive_use_any) %>%
    tidyr::pivot_longer(cols = -seq(1, ncolumns -3), names_to = "indicator", values_to = "value")
  return(samples_all)
}
