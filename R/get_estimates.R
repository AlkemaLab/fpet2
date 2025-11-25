
#' Get estimates for all-women model fits
#'
#' This function processes posterior samples from FPET to get estimates for
#' different marital status groups and indicators.
#' It can also handle subnational estimation if specified.
#'
#' @param samples tibble of posterior samples with columns for indicators, geo_unit_index, time_index, and value
#' @param geo_unit_index geo_unit_index from model fit
#' @param time_index time_index from model fit
#' @param subnational boolean, if TRUE, also computes national aggregates for subnational models
#' @param return_samples boolean, if TRUE, returns posterior samples along with summaries
#'
#' @returns tibble with estimates (mean and percentiles) for each indicator and marital status group
#' @keywords internal
get_estimates_fpem <- function(samples,
                               geo_unit_index,
                               time_index,
                               subnational = FALSE,
                               return_samples = FALSE){
  indicators_base <- c("unmet", "modern", "trad")
  groups <- c("", "unmarried_", "allwomen_")
  summ <- NULL
  samples_all <- list()
  for (i in 1:length(groups)){
    indicators <- paste0(groups[i], indicators_base)
    samples_some <- get_samples_for_indicators(indicators = indicators,
                                               samples = samples,
                                               geo_unit_index = geo_unit_index,
                                               time_index = time_index,
                                               stringtoremove = groups[i],
                                               rename_modern = TRUE)

    samples_all[[i]] <- summarize_fp_helper(samples_some) %>%
      dplyr::select(-c, -t) %>%
      mutate(marital_status =   case_when(
        groups[i] == "" ~ "married",
        groups[i] == "unmarried_" ~ "unmarried",
        groups[i] == "allwomen_" ~ "all")
      )
    summ <- bind_rows(
      summ,
      summarize_samples(samples_all[[i]]) %>%
        mutate(measure = "proportion"))
  }
  # if subnat, add aggregates
  if (subnational){ #"region_code" %in% geo_unit_index)
    indicators_base <- paste0(indicators_base, "_aggr")
    for (i in 1:3){
      list_counter <- length(groups) + i
      indicators <- paste0(groups[i], indicators_base)
      samples_some <- get_samples_for_indicators(indicators = indicators,
                                                 has_geo = FALSE,
                                                 samples = samples,
                                                 geo_unit_index = geo_unit_index,
                                                 time_index = time_index,
                                                 stringtoremove = groups[i],
                                                 rename_modern = TRUE) %>%
        # just remove groups[i] from indicators
        mutate(indicator = gsub(groups[i], "", indicator))
      samples_all[[list_counter]] <- summarize_fp_helper(samples_some) %>%
        dplyr::select( -t) %>%
        mutate(marital_status =   case_when(
          groups[i] == "" ~ "married",
          groups[i] == "unmarried_" ~ "unmarried",
          groups[i] == "allwomen_" ~ "all")
        ) %>%
        mutate(measure = "proportion", region_code = "National")
      summ <- bind_rows(
        summ,
        summarize_samples(samples_all[[list_counter]]) %>%
          mutate(measure = "proportion"))
    }}
  if (return_samples){
    return(list(
      posterior_samples = bind_rows(samples_all),
      estimates = summ
    ))
  } else {
    return(summ)
  }
}

#' Summarize samples
#'
#' @param samples tibble with `value` column and other grouping columns
#' @param percentiles vector of percentiles to compute
#'
#' @returns tibble with summary columns (mean and percentiles) for each group
#' @keywords internal
summarize_samples <- function(samples, # tibble with value column and lots of other info columns
                              percentiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)){
  group_cols <- names(samples)[names(samples) !="value"]
  summ <- samples %>%
    group_by(across(all_of(group_cols))) %>%
    dplyr::reframe(tibble::enframe(c(mean = mean(value),
                                     quantile(value, percentiles)), "percentile", "value"))
  return(summ)
}


#' Get modern from demand and demand satisfied
#'
#' Model returns Eta (demandsat) and d_Eta (demand)
#' # real demand = inv_tr_eta(tr_d_Eta_obs[geo_unit[i], time[i]]);
#' # modern[i] = inv_tr_eta(tr_Eta_obs[geo_unit[i], time[i]]) * demand;
#' # logit_unmetovernonmodern[i] = logit((demand - modern[i])/(1-modern[i]));

#' @param demand demand
#' @param demand_satisfied demand satisfied (modern/demand)
#'
#' @returns modern
#' @keywords internal
get_modern_fromdemandanddemandsat <- function(demand, demand_satisfied){
  demand * demand_satisfied
}

#' Get unmet from demand and modern
#'
#' @param demand demand
#' @param modern modern
#'
#' @returns unmet
#' @keywords internal
get_unmet_fromdemandandmodern <- function(demand, modern){
  demand - modern
}

#' Get estimates from global runs
#'
#'
#' Let's keep this function as we use it for the global runs
#' it worked for all women but not yet aggregates
#'
#' @param samples tibble of posterior samples with columns for indicators, geo_unit_index, time_index, and value
#' @param geo_unit_index geo_unit_index from model fit
#' @param time_index time_index from model fit
#' @param marital_status string, "married" or "unmarried"
#' @param return_samples boolean, if TRUE, returns posterior samples along with summaries
#' @param add_trad boolean, if TRUE, includes traditional use in the output
#'
#' @returns tibble with estimates (mean and percentiles) for each indicator and marital status group
#' @export
#'
get_estimates_globalruns <- function(samples, geo_unit_index,
                          time_index,
                          marital_status = "married", #"married", "unmarried",

                          # this is added as a column
                          return_samples = FALSE,
                          add_trad = TRUE){
  geo_unit <- ifelse(is.element("region_code", names(geo_unit_index)), "region_code", "iso")
  # iso needs to be updated to geo_unit
  group_cols <- c("c", "t", unique(c("iso", geo_unit)), "year", "indicator", "marital_status")
  if (!add_trad){ # not to use for all women
    indicators <- c("Eta", "d_Eta")
    samples_all <- get_samples_for_indicators(indicators = indicators, has_geo = TRUE,
                                              samples =samples,
                                              geo_unit_index =  geo_unit_index,
                                              time_index = time_index)
    samples_all <- samples_all %>%
      mutate(indicator = case_when(
        indicator == "d_Eta" ~ "demand",
        indicator == "Eta" ~ "demand_satisfied_modern",
        TRUE ~ indicator)) %>%
      tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
      mutate(
        modern = get_modern_fromdemandanddemandsat(demand = demand, demand_satisfied = demand_satisfied_modern),
        unmet = get_unmet_fromdemandandmodern(demand = demand, modern = modern)
      ) %>%
      tidyr::pivot_longer(cols = c(modern, unmet, demand_satisfied_modern, demand), names_to = "indicator", values_to = "value") %>%
      mutate(indicator = case_when(
        indicator == "modern" ~ "contraceptive_use_modern",
        indicator == "unmet" ~ "unmet_need_modern",
        #indicator == "trad" ~ "contraceptive_use_traditional",
        TRUE ~ indicator)) %>%
      mutate(marital_status = marital_status)
  } else {
    # if (marital_status == "unmarried" & all_women){
    #   indicators <- c("unmarried_d_Eta", "unmarried_unmet", "unmarried_trad")
    # } else {
      indicators <- c("d_Eta", "unmet", "trad")
   # }
    samples_all <- get_samples_for_indicators(indicators = indicators,
                                              samples = samples,
                                              has_geo = TRUE,
                                              geo_unit_index = geo_unit_index,
                                              time_index = time_index)
    samples_all <- samples_all  %>%
      mutate(indicator = case_when(
        indicator == "d_Eta" ~ "demand",
        TRUE ~ indicator)) %>%
      tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
      mutate(
        modern = demand - unmet,
        demand_satisfied_modern = modern/demand
      ) %>%
      tidyr::pivot_longer(cols = c(modern, unmet, demand_satisfied_modern, demand, trad), names_to = "indicator", values_to = "value") %>%
      mutate(indicator = case_when(
        indicator == "modern" ~ "contraceptive_use_modern",
        indicator == "unmet" ~ "unmet_need_modern",
        indicator == "trad" ~ "contraceptive_use_traditional",
        TRUE ~ indicator)) %>%
      mutate(marital_status = marital_status)
  }
  if (return_samples)
    return(samples_all)
  samples_all %>%
    group_by(across(all_of(group_cols))) %>%
    dplyr::reframe(tibble::enframe(c(mean = mean(value), quantile(value, c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))), "percentile", "value")) %>%
  mutate(measure = "proportion")
}

# summarize_fp_sample <- function(samples, # has demand, demand_satisfied, unmet, modern
#                                 #subnational,
#                                 #marital  = FALSE
# ) {
#   group_cols <- c("c", "t", "name_country", "year", "indicator")
#   #if (marital) {
#   #  group_cols <- c(group_cols, "marital_status")
#   #}
#   # # memory issues
#   # if (subnational) {
#   #   group_cols <- c(group_cols, "region_code")
#   # }
#   # probaly a more efficient solution but this works
#   summarize_fp_helper <- function(samples){
#     ncolumns <- dim(samples)[2]
#     estimates <- samples %>%
#       tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
#       dplyr::select(-draw) %>%
#       mutate(contraceptive_use_any = contraceptive_use_modern + contraceptive_use_traditional,
#              non_use = 1 - contraceptive_use_any,
#              unmet_need_modern = unmet_need_any + contraceptive_use_traditional,
#              demand = contraceptive_use_any + unmet_need_any,
#              # demand_modern added here for consistency with fpet1, note that it is just demand so could be removed
#              demand_modern = demand,
#              demand_satisfied = contraceptive_use_any/demand_modern,
#              demand_satisfied_modern = contraceptive_use_modern/demand_modern,
#              no_need = 1 - demand,
#              ratio_modern_any = contraceptive_use_modern/contraceptive_use_any) %>%
#       tidyr::pivot_longer(cols = -seq(1, ncolumns -3), names_to = "indicator", values_to = "value") %>%
#       group_by(across(all_of(group_cols))) %>%
#       dplyr::reframe(tibble::enframe(c(mean = mean(value), quantile(value, c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))), "percentile", "value"))
#     return(estimates)
#   }
#   if (!subnational){
#     estimates <-  summarize_fp_helper(samples)
#   } else {
#     estimates <-  NULL
#     for (region in unique(samples$region_code)){
#       estimates <- bind_rows(estimates,
#                              summarize_fp_helper(samples %>% dplyr::filter(region_code == region)) %>% mutate(region_code = region))
#     }
#   }
#   # group_by(across(all_of(group_cols))) %>%
#   # group_modify(~ {
#   # add mean
#   #   quantile(.x$value, probs =  c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)) %>%
#   #     tibble::enframe(name = "percentile", value = "value")
#   #   })
#   return(estimates)
# }

# get_samples_fpindicators <- function(fit_data, samples, subnational, all_women = FALSE) {
#   # unmet here refers to unmet_modern
#   results <- get_samples_fpindicators_one_group(fit_data, samples, agg = FALSE, all_women = all_women)
#
#   if (subnational) {
#     results <- bind_rows(
#       results,
#       get_samples_fpindicators_one_group(fit_data, samples, agg = TRUE, all_women = all_women) %>%
#         mutate(c = 0)
#     )
#   }
#
#   return(results)
# }
#
# get_samples_fpindicators_one_group <- function(fit_data, samples,  agg, all_women) {
#
#   # unmet here refers to unmet_modern
#   indicators <- c("modern", "unmet", "trad")
#   if (agg) {
#     indicators <- paste0(indicators, "_agg")
#     geo_unit_index <- fit_data$agg_geo_unit_index |>
#       rename(
#         name_country = agg_unit_name,
#         c = i) |>
#       mutate(
#         region_code = "National"
#       )
#   } else {
#     geo_unit_index <- fit_data$geo_unit_index
#   }
#
#
#   if (!all_women){
#     parameter_pattern <- paste0(
#       "^(",
#       paste(indicators, collapse = "|"),
#       ")\\[(\\d+),(\\d+)\\]$")
#
#     samples$draws(indicators) |>
#       posterior::as_draws_df() |>
#       tidyr::pivot_longer(cols = !(.chain:.draw)) |>
#       dplyr::rename_with(~ stringr::str_replace(.x, "^\\.", "")) |>
#       dplyr::mutate(
#         stringr::str_match(name, parameter_pattern) |>
#           as_tibble(.name_repair = ~ c("matched", "parameter", "c", "t"))
#       ) |>
#       dplyr::mutate(
#         c = as.integer(c),
#         t = as.integer(t)
#       ) |>
#       dplyr::inner_join(geo_unit_index %>% select(- any_of(c("is_unmarried_sexual_activity", "hier_regional_unmarried",
#                                                              "cluster", "subcluster"))), by = "c") |>
#       dplyr::inner_join(fit_data$time_index, by = "t") %>%
#       dplyr::select(-chain, -iteration,  - name, -matched) %>%
#       dplyr::rename(indicator = parameter) %>%
#       dplyr::mutate(indicator = gsub("_agg", "", indicator)) %>%
#       dplyr::tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
#       dplyr::mutate(unmet_need_any = unmet - trad) %>%
#       dplyr::select(-unmet) %>%
#       tidyr::pivot_longer(cols = c(modern, trad, unmet_need_any), names_to = "indicator", values_to = "value") %>%
#       dplyr::mutate(indicator = case_when(
#         indicator == "modern" ~ "contraceptive_use_modern",
#         indicator == "trad" ~ "contraceptive_use_traditional",
#         TRUE ~ indicator))
#   } else {
#     # for all women
#     # copy-pasting for starters to not break the first option
#
#     # add extra dimension to get married and unmarried
#     parameter_pattern <- paste0(
#       "^(",
#       paste(indicators, collapse = "|"),
#       ")\\[(\\d+),(\\d+),(\\d+)\\]$")
#
#     #fit_all$samples$draws(indicators) |>
#     samples$draws(indicators) |>
#       posterior::as_draws_df() |>
#       tidyr::pivot_longer(cols = !(.chain:.draw)) |>
#       dplyr::rename_with(~ stringr::str_replace(.x, "^\\.", "")) |>
#       dplyr::mutate(
#         stringr::str_match(name, parameter_pattern) |>
#           as_tibble(.name_repair = ~ c("matched", "parameter", "m", "c", "t"))
#       ) |>
#       dplyr::mutate(
#         m = as.integer(m),
#         c = as.integer(c),
#         t = as.integer(t)
#       ) |>
#       dplyr::inner_join(geo_unit_index %>% select(- any_of(c("is_unmarried_sexual_activity", "hier_regional_unmarried",
#                                                              "cluster", "subcluster"))), by = "c") |>
#       dplyr::inner_join(fit_data$time_index, by = "t") %>%
#       dplyr::select(-chain, -iteration,  - name, -matched) %>%
#       dplyr::rename(indicator = parameter) %>%
#       dplyr::mutate(indicator = gsub("_agg", "", indicator)) %>%
#       tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
#       dplyr::mutate(unmet_need_any = unmet - trad) %>%
#       dplyr::select(-unmet) %>%
#       tidyr::pivot_longer(cols = c(modern, trad, unmet_need_any), names_to = "indicator", values_to = "value") %>%
#       dplyr::mutate(indicator = case_when(
#         indicator == "modern" ~ "contraceptive_use_modern",
#         indicator == "trad" ~ "contraceptive_use_traditional",
#         TRUE ~ indicator))
#
#
#   }
#
# }
#
#
#
#
