
#' Get residuals
#'
#' @param fit fit object from a global fit
#' @param is_married boolean
#'
#' @returns tibble with survey data and residuals
#' @keywords internal
get_residuals <- function(fit, is_married){
  # some additional stuff outputted
  estimates_samples <- get_estimates_globalruns(samples = fit$samples,
                                     geo_unit_index = fit$geo_unit_index,
                                     time_index = fit$time_index,
                                     marital_status = ifelse(is_married, "married", "unmarried"),
                                     return_samples = TRUE,
                                     add_trad =  ifelse(fit$runstep %in% c("step1a", "step1b"), FALSE, TRUE ))
  est_samples_use <-
    estimates_samples %>%
    dplyr::mutate(indicator = paste0("est_", indicator)) %>%
    tidyr::pivot_wider(names_from = indicator, values_from = value)
  est_samples_use

  # watch out: 100s imputed!!!
  dat_use <- tibble(
    i = seq(1, fit$stan_data$N),
    DM1_y = ifelse(fit$stan_data$DM1_obs_isna, NA, fit$stan_data$DM1_y),
    DM2_y = ifelse(fit$stan_data$DM2_obs_isna, NA, fit$stan_data$DM2_y),
    modern = inv_logit(DM1_y),
    unmet = inv_logit(DM2_y)*(1-modern),
    demand = modern + unmet,
    demand_satisfied_modern = modern/demand)
  dat_use

  geo_unit_index <- fit$geo_unit_index
  time_index <- fit$time_index
  samples <- fit$samples
  indicators <- c("DM1_scale", "DM2_scale")
  parameter_pattern <- paste0(
    "^(",
    paste(indicators, collapse = "|"),
    ")\\[(\\d+)\\]$")
  samples_res <-
    samples$draws(indicators) |>
    posterior::as_draws_df() |>
    tidyr::pivot_longer(cols = !(.chain:.draw)) |>
    # dplyr::rename_with(~ stringr::str_replace(.x, "^\\.", "")) |>
    dplyr::mutate(
      stringr::str_match(name, parameter_pattern) |>
        as_tibble(.name_repair = ~ c("matched", "parameter", "i"))) |>
    dplyr::mutate(
      i = as.integer(i)
    ) |>
    dplyr::rename(indicator = parameter) %>%
    dplyr::select(-.chain, -.iteration, - name, -matched) %>%
    dplyr::left_join( tibble(
      i = seq(1, fit$stan_data$N),
      c = fit$stan_data$geo_unit, t = fit$stan_data$time), by = "i") %>%
    tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
    dplyr::left_join(est_samples_use, by = c("c", "t", ".draw" = "draw")) %>%
    dplyr::left_join(dat_use)
  summ_res <- samples_res  %>%
    dplyr::mutate(st_res_modern = (DM1_y - logit(est_contraceptive_use_modern))/DM1_scale,
           st_res_unmetovernonmodern = (DM2_y -
                                          logit(est_unmet_need_modern/(1-est_contraceptive_use_modern)))/DM2_scale,
           res_demand = demand - est_demand,
           res_demandsat = demand_satisfied_modern - est_demand_satisfied_modern,
           res_modern = modern - est_contraceptive_use_modern)
  summ_res2 <-
    summ_res %>%
    dplyr::group_by(i) %>%
    dplyr::summarise(
      mean_st_res_modern = mean(st_res_modern, na.rm = TRUE),
      mean_st_res_unmetovernonmodern = mean(st_res_unmetovernonmodern),
      mean_res_demand = mean(res_demand),
      mean_res_demandsat = mean(res_demandsat),
      mean_res_modern = mean(res_modern),
      # to keep to use
      mean_modern = mean(est_contraceptive_use_modern),
      mean_unmet = mean(est_unmet_need_modern),
      mean_demand = mean(est_demand),
      mean_demand_satisfied_modern = mean(est_demand_satisfied_modern),
      c = first(c), t = first(t)
    ) %>%
    dplyr::left_join(fit$geo_unit_index)
  return(summ_res2)
}
