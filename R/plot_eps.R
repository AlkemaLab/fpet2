
#' Plot Epsilon Innovations
#'
#' @param fit fit object from a global fit
#' @param is_married boolean
#' @param estimates_df dataframe of estimates
#' @param add_trad boolean, whether to add traditional epsilon innovations
#'
#' @returns NULL, saves plots to fit$output_dir/eps.pdf
#' @keywords internal
plot_eps <- function(fit, is_married, estimates_df,
                     add_trad = FALSE){



  # estimates_df <- get_estimates(samples = fit$samples,
  #                               geo_unit_index = fit$geo_unit_index,
  #                               time_index = fit$time_index,
  #                               is_married = is_married)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot. Please install it.", call. = FALSE)
  }
  geo_unit_index <- fit$geo_unit_index
  time_index <- fit$time_index
  samples <- fit$samples
  if (add_trad){
    indicators <- c("Epsilon_innovation", "d_Epsilon_innovation", "z_Epsilon_innovation")
  } else {
    indicators <- c("Epsilon_innovation", "d_Epsilon_innovation")
  }
  # iso needs to be updated to geo_unit
  group_cols <- c("c", "t", "iso", "year", "indicator")
  parameter_pattern <- paste0(
    "^(",
    paste(indicators, collapse = "|"),
    ")\\[(\\d+),(\\d+)\\]$")
  res <-
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
    rename(indicator = parameter) %>%
    # mutate(indicator = case_when(
    #   indicator == "d_Epsilon_innovation" ~ "eps_inno_demand",
    #   indicator == "Epsilon_innovation" ~ "demand_satisfied_modern",
    #   TRUE ~ indicator)) %>%
    group_by(c,t, indicator) %>%
    summarise(eps_innov = mean(value),
              abs_eps_innov = abs(eps_innov)) %>%
    dplyr::inner_join(geo_unit_index #%>% select(- any_of(c("is_unmarried_sexual_activity", "hier_regional_unmarried",
                      #                     "cluster", "subcluster")))
                      , by = "c") |>
    dplyr::inner_join(time_index, by = "t") %>%
    ungroup() %>%
    left_join(tibble(c = seq(1, nrow(geo_unit_index)),
                     t_min = c(fit$stan_data$t_min),
                     t_max = c(fit$stan_data$t_max))) %>%
    dplyr::filter(t >= t_min & t <= t_max) %>%
    left_join(
      estimates_df %>%
        dplyr::filter(percentile == "mean") %>%
        tidyr::pivot_wider(names_from = indicator, values_from = value) )
  #names(res)
  if (!is_married){
    res <- res %>%
      rename(facet_region = is_unmarried_sexual_activity)
  } else {
    res <- res %>%
      rename(facet_region = cluster)
  }
  #return(res)

  pdf(file = file.path(fit$output_dir, "eps.pdf"), width = 12, height = 8)
  for (eps in unique(res$indicator)){
    for (xvar in c("demand_satisfied_modern", "demand", "year", "contraceptive_use_modern")){
      p <- res %>%
        dplyr::filter(indicator == eps) %>%
        ggplot(ggplot2::aes(x = !!sym(xvar), y = eps_innov, color = facet_region)) +
        geom_point() +
        geom_smooth(method = "loess", color = "red", se = TRUE,
                    method.args = list(family="symmetric", span = 0.5)) +
        ggtitle(paste(eps)) +
        theme(legend.position = "none")
    print(p)

      p <- res %>%
        dplyr::filter(indicator == eps) %>%
        ggplot(ggplot2::aes(x = !!sym(xvar), y = eps_innov)) +
        geom_point() +
        geom_smooth(method = "loess", color = "red", se = TRUE,
                    method.args = list(family="symmetric", span = 0.5)) +
        ggtitle(paste(eps)) +
        theme(legend.position = "none")+
        facet_wrap(~ facet_region, scales = "free" )
      print(p)
    }
  }

  dev.off()
  return(NULL)
}
