

#' Plot residuals
#'
#' @param fit fit object from a global fit
#' @param is_married boolean
#'
#' @returns NULL, but saves residual plots to output_dir
#' @keywords internal
plot_residuals <- function(fit, is_married){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot. Please install it.", call. = FALSE)
  }

  # for fit1a, get residuals based on observed demand and demand satisfied..
  # note that these are only a subset of the observations
  # also get them for modern and unmet

  # // get means for all, to use for pma and non-pma
  # vector[N] modern;
  # vector[N] logit_unmetovernonmodern;
  # for (i in 1:N){
  #   real demand = inv_tr_eta(tr_d_Eta_obs[geo_unit[i], time[i]]);
  #   modern[i] = inv_tr_eta(tr_Eta_obs[geo_unit[i], time[i]]) * demand;
  #   logit_unmetovernonmodern[i] = logit((demand - modern[i])/(1-modern[i]));
  # }
  #
  # // non-pma data
  # for(i in 1:N) {
  #   if(held_out[i] == 0 && ispma[i] == 0) {
  #     if(DM1_obs_isna[i] == 0) {
  #       DM1_y[i] ~ normal(logit(modern[i]), DM1_scale[i]);
  #     }
  #     if(DM2_obs_isna[i] == 0) {
  #       DM2_y[i] ~ normal(logit_unmetovernonmodern[i], DM2_scale[i]);
  #     }
  #   }
  #devtools::load_all(here::here())

  summ_res2 <- get_residuals(fit, is_married)

  if (!is_married){
    summ_res2 <- summ_res2 %>%
      rename(facet_region = is_unmarried_sexual_activity)
  } else {
    summ_res2 <- summ_res2 %>%
      rename(facet_region = cluster)
  }


  pdf(file = file.path(fit$output_dir, "residuals.pdf"), width = 12, height = 8)

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = mean_modern, y = mean_st_res_modern, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = mean_unmet, y = mean_st_res_unmetovernonmodern, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = mean_demand, y = mean_res_demand, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = mean_demand_satisfied_modern, y = mean_res_demandsat, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = mean_modern, y = mean_res_modern, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = t, y = mean_res_modern, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = t, y = mean_res_demandsat, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))

  print(summ_res2 %>%
          ggplot2::ggplot(ggplot2::aes(x = t, y = mean_res_demand, color = facet_region)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "loess", color = "red", se = TRUE,
                method.args = list(family="symmetric", span = 0.5)
    ))
  dev.off()
  return(NULL)
}
