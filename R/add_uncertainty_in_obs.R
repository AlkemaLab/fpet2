
#' Calculate uncertainty in observations for all-women fit
#'
#' For all-women FPET fit, adds uncertainty in observations based on sampling errors and, if specified, also adds NSE uncertainty.

#' @param samples posterior::draws object from cmdstanr fit
#' @param data_allwomen list with data for all women used in model fitting
#' @param is_married boolean
#' @param perc_low percentile for lower bound of uncertainty interval, defaults to 0.025
#' @param perc_up percentile for upper bound of uncertainty interval, defaults to 0.975
#' @param add_nse_uncertainty boolean, if TRUE (default), adds NSE uncertainty in addition to sampling error uncertainty
#'
#' @returns survey data for both marital groups with columns for point estimates and uncertainty intervals added
#' @keywords internal
add_uncertainty_in_obs_allwomen <- function(samples,
                                       data_allwomen,
                                       is_married = TRUE,
                                       perc_low = 0.025, perc_up = 0.975,
                                       add_nse_uncertainty = TRUE){

  tmp <- list()
  marital_status <- ifelse(is_married, "married", "unmarried")
  tmp$data <- data_allwomen[[marital_status]]$data
  tmp$stan_data <- data_allwomen[[marital_status]]$stan_data # uses DMx_y, DMx_s, DMnotobs
  if (is_married){
    samples_scale <- list()
    samples_scale[["DM1_scale"]] <- samples$draws("DM1_scale")
    samples_scale[["DM2_scale"]] <- samples$draws("DM2_scale")
    samples_scale[["DM3_scale"]] <- samples$draws("DM3_scale")
    res <- add_uncertainty_in_obs(tmp,
                                  samples_scale = samples_scale,
                                  perc_low, perc_up,
                                  add_nse_uncertainty)
  } else {
    # assuming we have trad
    samples_scale <- list()
    samples_scale[["DM1_scale"]] <- samples$draws("unmarried_DM1_scale")
    samples_scale[["DM2_scale"]] <- samples$draws("unmarried_DM2_scale")
    samples_scale[["DM3_scale"]] <- samples$draws("unmarried_DM3_scale")
    res <- add_uncertainty_in_obs(tmp,
                                  samples_scale = samples_scale,
                                  perc_low, perc_up,
                                  add_nse_uncertainty)
  }
  return(res)
}


#' Calculate uncertainty in observations for fit to one marital group
#'
#' For one marital group, adds uncertainty in observations based on sampling errors and, if specified, also adds NSE uncertainty.
#' To be used for global fits for one marital group.


#' @param fit list with fit information, including data and posterior samples
#' @param samples_scale if provided, used instead of fit$samples$draws("DMX_scale")
#' @param is_married boolean
#' @param perc_low percentile for lower bound of uncertainty interval, defaults to 0.025
#' @param perc_up percentile for upper bound of uncertainty interval, defaults to 0.975
#' @param add_nse_uncertainty boolean, if TRUE (default), adds NSE uncertainty in addition to sampling error uncertainty
#'
#' @returns survey data with columns for point estimates and uncertainty intervals added
#' @export
add_uncertainty_in_obs <- function(fit, # to get data and sampling errors
                                   samples_scale = NULL, # if provided, used instead of fit$samples$draws("DMX_scale")
                                   # list with samples_scale[["DM1_scale"]] etc
                                   perc_low = 0.025, perc_up = 0.975,
                                   add_nse_uncertainty = TRUE){
  dat = fit$data
  DM1_y = fit$stan_data$DM1_y
  DM1_s = fit$stan_data$DM1_s
  DM2_y = fit$stan_data$DM2_y
  DM2_s = fit$stan_data$DM2_s
  if (!is.null(fit$stan_data$DM3_y)){
    add_trad = TRUE
    DM3_y = fit$stan_data$DM3_y
    DM3_s = fit$stan_data$DM3_s
  } else {
    add_trad = FALSE
  }
  nobs <- dim(dat)[1]
  if (nobs != length(DM1_y)){
    stop("nobs in dat and y differ")
  }

  if (add_nse_uncertainty){
    if (is.null(samples_scale)){
      samples <- fit$samples
      samples_scale <- list()
      samples_scale[["DM1_scale"]] <- samples$draws("DM1_scale")
      samples_scale[["DM2_scale"]] <- samples$draws("DM2_scale")
      if (add_trad){
        samples_scale[["DM3_scale"]] <- samples$draws("DM3_scale")
      }
    }
    DM1_scale = t(get_scales_si(scale_samples = samples_scale[["DM1_scale"]]))
    DM2_scale = t(get_scales_si(scale_samples = samples_scale[["DM2_scale"]]))
    if (add_trad){
      DM3_scale = t(get_scales_si(scale_samples = samples_scale[["DM3_scale"]]))
    }
  } else {
    print("We only assess the uncertainty based on sampling errors")
    # just make into matrix such that format is the same
    DM1_scale <- matrix(DM1_s, nobs, 1)
    if (nobs != dim(DM1_scale)[1]){
      stop("nobs in dat and scale_y differ")
    }
    DM2_scale <- matrix(DM2_s, nobs, 1)
    if (add_trad){
      DM3_scale = matrix(DM2_s, nobs, 1)
    }
  }
  nsamples_scale <- dim(DM1_scale)[2]
  nsamples <- ifelse(nsamples_scale == 1, 2000, nsamples_scale)
  set.seed(12345)
  # simulate observations
  #print("internal note: to be extended if we add obs on total use only")
  mcpr_is <-   unmet_need_modern_is <- trad_is <- matrix(NA, nobs, nsamples)

  # for dm3
  # data <- data %>%
  #   mutate(logit_trad_overunmet = ifelse(DM2_obs_isna, logit_trad_overnotmodern, logit_trad_overunmet),
  #          se_logit_trad_overunmet = ifelse(DM2_obs_isna, se_logit_trad_overnotmodern, se_logit_trad_overunmet))


  for (i in 1:nobs){
    trad_is[i,] <- dat$contraceptive_use_traditional[i]
    if (fit$stan_data$DM1_obs_isna[i]==1){
      # other obs not used either.. so reset to NA
      mcpr_is[i,] <- NA
      trad_is[i,] <- NA
      unmet_need_modern_is[i,] <- NA
    } else {
      mcpr_is[i,] <- inv_logit(rnorm(nsamples, DM1_y[i], DM1_scale[i,]))
      if (fit$stan_data$DM2_obs_isna[i]==0){ # if we have unmet_modern
        unmet_need_modern_is[i,] <- inv_logit(rnorm(nsamples, DM2_y[i], DM2_scale[i,]))*(1-median(mcpr_is[i,]))
        if (add_trad){
          # note that adding sampled uncertainty is not accurate here either
          if (fit$stan_data$DM3_obs_isna[i]==0){
            trad_is[i,] <- inv_logit(rnorm(nsamples, DM3_y[i], DM3_scale[i,]))*median(unmet_need_modern_is[i,])
          }
        }
      } else {
        unmet_need_modern_is[i,] <- NA
        if (add_trad){
          if (fit$stan_data$DM3_obs_isna[i]==0){
            trad_is[i,] <- inv_logit(rnorm(nsamples, DM3_y[i], DM3_scale[i,]))*(1-median(mcpr_is[i,]))
          }
        }
      } # end if DM2_obs_isna
    } # end mcpr not missing
  }

  any_is <- mcpr_is + trad_is
  unmet_need_any_is <- unmet_need_modern_is - trad_is
  demand_is <- mcpr_is + unmet_need_modern_is
  demand_satisfied_modern_is <- mcpr_is/demand_is
  demand_satisfied_any_is <- any_is/demand_is
  ratio_modern_total_is <- mcpr_is/any_is

  dat %>%
    dplyr::mutate(
      low_contraceptive_use_modern = apply(mcpr_is, 1, quantile, perc_low, na.rm = TRUE),
      up_contraceptive_use_modern = apply(mcpr_is, 1, quantile, perc_up, na.rm = TRUE),
      est_contraceptive_use_modern = contraceptive_use_modern,
      low_contraceptive_use_traditional = apply(trad_is, 1, quantile, perc_low, na.rm = TRUE),
      up_contraceptive_use_traditional = apply(trad_is, 1, quantile, perc_up, na.rm = TRUE),
      est_contraceptive_use_traditional = contraceptive_use_traditional,
      low_contraceptive_use_any = apply(any_is, 1, quantile, perc_low, na.rm = TRUE),
      up_contraceptive_use_any = apply(any_is, 1, quantile, perc_up, na.rm = TRUE),
      est_contraceptive_use_any = contraceptive_use_any,
      low_unmet_need_any = apply(unmet_need_any_is, 1, quantile, perc_low, na.rm = TRUE),
      up_unmet_need_any = apply(unmet_need_any_is, 1, quantile, perc_up, na.rm = TRUE),
      est_unmet_need_any = unmet_need_any,
      low_demand = apply(demand_is, 1, quantile, perc_low, na.rm = TRUE),
      up_demand = apply(demand_is, 1, quantile, perc_up, na.rm = TRUE),
      est_demand = demand,
      low_demand_satisfied_modern = apply(demand_satisfied_modern_is, 1, quantile, perc_low, na.rm = TRUE),
      up_demand_satisfied_modern = apply(demand_satisfied_modern_is, 1, quantile, perc_up, na.rm = TRUE),
      est_demand_satisfied_modern = demand_satisfied_modern,
      low_demand_satisfied_any = apply(demand_satisfied_any_is, 1, quantile, perc_low, na.rm = TRUE),
      up_demand_satisfied_any = apply(demand_satisfied_any_is, 1, quantile, perc_up, na.rm = TRUE),
      est_demand_satisfied_any = demand_satisfied_any,
      low_ratio_modern_any = apply(ratio_modern_total_is, 1, quantile, perc_low, na.rm = TRUE),
      up_ratio_modern_any = apply(ratio_modern_total_is, 1, quantile, perc_up, na.rm = TRUE),
      est_ratio_modern_any = ratio_modern_any,
      low_unmet_need_modern = apply(unmet_need_modern_is, 1, quantile, perc_low, na.rm = TRUE),
      up_unmet_need_modern = apply(unmet_need_modern_is, 1, quantile, perc_up, na.rm = TRUE),
      est_unmet_need_modern = unmet, # unmet is unmet_modern in data
      # added for validation
      held_out = fit$stan_data$held_out
    )
}

get_scales_si <- function(scale_samples){
  parameter_pattern <- "^(DM1_scale|DM2_scale|DM3_scale|unmarried_DM1_scale|unmarried_DM2_scale|unmarried_DM3_scale)\\[(\\d+)\\]$"
  scale_samples |>
    posterior::as_draws_df() |>
    tidyr::pivot_longer(cols = !(.chain:.draw)) |>
    dplyr::rename_with(~ stringr::str_replace(.x, "^\\.", "")) |>
    dplyr::mutate(
      stringr::str_match(name, parameter_pattern) |>
        as_tibble(.name_repair = ~ c("matched", "parameter", "i"))
    ) |>
    dplyr::select(draw, value, i) %>%
    tidyr::pivot_wider(names_from = i) %>%
    dplyr::select(-draw)
}
