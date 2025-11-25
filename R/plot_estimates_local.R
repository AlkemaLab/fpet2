#' Plot FP estimates
#'
#' Plot FP estimates for multiple countries/regions, marital status, and indicators
#'
#' @param iso_select name or names of country iso or region code
#' @param subnational logical, whether subnational or national
#' @param results output from fit_fpem()
#' @param marital_status_all choice out of c("married", "unmarried", "all")
#' @param indicator_select_all choice out of c("contraceptive_use_modern", "unmet_need_modern",
#' "demand", "demand_satisfied_modern")
#' @param dat_emu data frame of emu data, if NULL, will not plot emu data
#' @param results2 optional second model results to compare
#' @param results3 optional third model results to compare
#' @param results4 optional fourth model results to compare
#' @param modelnames names of models to use in legend if comparing more than one model
#' @param add_title logical, whether to add a title to the plot
#' @param nrow_plot when plotting 4 indicators, use 1 or 2 rows
#' @param save_plots logical, whether to save plots as pdf
#' @param cols_emus colors to use for emu data types
#' @param output_folder folder to save plots in
#' @param plot_name name of pdf file to save plots in
#'
#' @returns list of plots for selected countries/regions, marital status, and indicators
#' @export
#'
plot_estimates_local_all <- function( iso_select = NULL,
                                      subnational = FALSE,
                                     results,
                                     marital_status_all =  c("married", "unmarried", "all"),
                                     indicator_select_all = c(
                                       "contraceptive_use_modern",
                                       "unmet_need_modern", "demand", "demand_satisfied_modern",
                                       "contraceptive_use_traditional",
                                       "contraceptive_use_any",
                                       "unmet_need_any",
                                       "demand_satisfied"),
                                     dat_emu = NULL,
                                     results2 = NULL,
                                     results3 = NULL,
                                     results4 = NULL,
                                     modelnames = NULL,
                                     add_title = TRUE,
                                     nrow_plot = 2,
                                     save_plots = TRUE,
                                     cols_emus = c("clients"    = "red",        "visits"   = "blue",
                                                   "facilities" = "darkgreen", "users"  = "orange"),
                                     output_folder = NULL,
                                     plot_name = "estimates"){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package 'ggnewscale' is required for this plot. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for this plot. Please install it.", call. = FALSE)
  }


  #if (is.null(plot_title)){
  #  plot_title <- std_country
  #}
  plot_title <- "Fit"
  all_plots <- list()
  geo_col <- ifelse(subnational, "region_code", "iso")


  # Get unique country or subnational codes
  if (is.null(iso_select)){
    model_country_codes <- results$estimates[[geo_col]] %>% unique()
  } else {
    model_country_codes <- iso_select
  }
  country_codes <- model_country_codes
  fit_data <- results$observations

  for (i in 1:length(country_codes)) {
    filtered_data_all <- fit_data %>%
      dplyr::filter(.data[[geo_col]] == country_codes[i])
    estimates_all <- results$estimates %>%
      dplyr::filter(.data[[geo_col]] == country_codes[i])

    estimates2_all <- if (!is.null(results2)) results2$estimates %>%
      dplyr::filter(.data[[geo_col]] == country_codes[i]) else NULL
    estimates3_all <- if (!is.null(results3)) results3$estimates %>%
      dplyr::filter(.data[[geo_col]] == country_codes[i]) else NULL
    estimates4_all <- if (!is.null(results4)) results4$estimates %>%
      dplyr::filter(.data[[geo_col]] == country_codes[i]) else NULL
    all_plots[[country_codes[i]]] <- list()

    for (marital_status_select in marital_status_all){
      pcount <- 1
      plots <- list()
      for (indicator_select in indicator_select_all){
        if (marital_status_select != "all" &
            indicator_select %in% names(results$observations)){
          filtered_data  <- filtered_data_all %>%
            dplyr::filter(
              is_in_union == ifelse(marital_status_select == "married", "Y", "N")) %>%
            dplyr::rename(indicator = paste0("est_", indicator_select),
                   indicator_low = paste0("low_", indicator_select),
                   indicator_up = paste0("up_", indicator_select)
            )


          if (all(is.na(filtered_data$indicator))){
            filtered_data  <- NULL
          } else {
            filtered_data  <- filtered_data %>%
              dplyr::filter(!is.na(indicator))
          }
        } else {
          filtered_data  <- NULL
        }
        estimates <- estimates_all %>%
          dplyr::filter(measure == "proportion",
                 indicator == indicator_select,
                 marital_status == marital_status_select)
        if (!is.null(results2)){
          estimates2 <- estimates2_all %>%
            dplyr::filter(measure == "proportion",
                   indicator == indicator_select,
                   marital_status == marital_status_select)
        } else {
          estimates2 <- NULL
        }
        if (!is.null(results3)){
          estimates3 <- estimates3_all %>%
            dplyr::filter(measure == "proportion",
                   indicator == indicator_select,
                   marital_status == marital_status_select)
        } else {
          estimates3 <- NULL
        }
        if (!is.null(results4)){
          estimates4 <- estimates4_all %>%
            dplyr::filter(measure == "proportion",
                   indicator == indicator_select,
                   marital_status == marital_status_select)
        } else {
          estimates4 <- NULL
        }
        plots[[pcount]] <- plot_estimates_local(estimates,
                                                filtered_data,
                                                indicator_select,
                                                estimates2 = estimates2,
                                                estimates3 = estimates3,
                                                estimates4 = estimates4,
                                                modelnames = modelnames)
        if (!is.null(results$dat_emu) & indicator_select == "contraceptive_use_modern"){
          dat_emu_select <-
            results$dat_emu %>%
            dplyr::filter(.data[[geo_col]] == country_codes[i]) %>%
            dplyr::filter(case_when(
                marital_status_select == "all" & emu_for_allwomen_j == 1 ~ TRUE,
                marital_status_select == "married" & emu_for_allwomen_j == 0 ~ TRUE,
                TRUE ~ FALSE))
          if (dim(dat_emu_select)[1] > 0){
            plots[[pcount]] <- plots[[pcount]] +
              ggnewscale::new_scale_color() +
              ggplot2::geom_errorbar(data = dat_emu_select,
                                    ggplot2::aes(x = year, y = emu, ymin = emu_lower, ymax = emu_upper, color = ss_type), alpha = 0.3) +
              ggplot2::geom_point(data = dat_emu_select,
                         ggplot2::aes(x = year, y = emu, color = ss_type)) +
              ggplot2::scale_colour_manual(values = cols_emus)


          }
        }
        #print( plots[[pcount]])
        pcount <- pcount + 1
      }
      if (length(indicator_select_all) !=1 ){
        if (length(indicator_select_all) == 4){
          if (nrow_plot == 2){
            p <- ggpubr::ggarrange(plotlist = plots,
                           ncol = 2, nrow = 2, common.legend = TRUE,
                           legend = "bottom")#, title = std_country)
          } else {
            p <- ggpubr::ggarrange(plotlist = plots,
                                   ncol = 4, nrow = 1, common.legend = TRUE,
                                   legend = "bottom")#, title = std_country)
          }
        } else {
          if (length(indicator_select_all) == 8){
            p <- ggpubr::ggarrange(plotlist = plots,
                           ncol = 4, nrow = 2, common.legend = TRUE,
                           legend = "bottom")#, title = std_country)
          } else {
            if (length(indicator_select_all) == 3){
              p <- ggpubr::ggarrange(plotlist = plots,
                             ncol = 3, nrow = 1, common.legend = TRUE,
                             legend = "bottom")#, title = std_country)
            } else {
              stop("implement this!")
        }}}

      } else {
        p <- ggpubr::ggarrange(plotlist = plots,
                       ncol = 1, nrow = 1, common.legend = TRUE,
                       legend = "right")#, title = std_country)
      }
      if (add_title){
        p <- ggpubr::annotate_figure(p,
                                     top = ggpubr::text_grob(paste0(country_codes[i], " ", marital_status_select),
                                                face = "bold", size = 14))
      }
      all_plots[[country_codes[i]]][[marital_status_select]] <- p
    }# end marital
  } # end countries
  if (save_plots) {
    if (!is.null(output_folder)) {
      output_dir <- output_folder
    #} else if (dir.exists(results$output_dir)) {
    #  output_dir <- results$output_dir
    } else {
      stop("Please provide a valid output_folder to save the plot in.")
    }
    pdf(file.path(output_dir, paste0(plot_name,".pdf")), width = 11, height = 6)
    for (plot in all_plots) {
      print(plot)
    }
    dev.off()
  }
  return(all_plots)
}


#' Plot FP estimates for one country/region, marital status, and indicator
#'
#' @param estimates estimates
#' @param filtered_data observations
#' @param indicator_select indicator
#' @param estimates2 optional second model estimates to compare
#' @param estimates3 optional third model estimates to compare
#' @param estimates4 optional fourth model estimates to compare
#' @param modelnames names of models to use in legend if comparing more than one model
#' @param cols_sourcetypes colors to use for data source types
#'
#' @returns ggplot object
#' @keywords internal
plot_estimates_local <- function(estimates, filtered_data,
                                 indicator_select,
                                 estimates2 = NULL,
                                 estimates3 = NULL,
                                 estimates4 = NULL,

                                 modelnames = c("model1", "model2"),
                                 cols_sourcetypes = c("DHS" = "red", "DHS0" = "red", "MICS" = "blue",
                                                      "PMA" = "darkgreen", "Other" = "orange", "National survey" = "purple")){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot. Please install it.", call. = FALSE)
  }
  model <- NULL
  if (!is.null(estimates2)){
    estimates <- bind_rows(
      estimates %>% mutate(model = modelnames[1]),
      estimates2 %>% mutate(model = modelnames[2])
    )
  }
  if (!is.null(estimates3)){
    estimates <- bind_rows(
      estimates,
      estimates3 %>% mutate(model = modelnames[3])
    )
  }
  if (!is.null(estimates4)){
    estimates <- bind_rows(
      estimates,
      estimates4 %>% mutate(model = modelnames[4])
    )
  }
  p <- estimates %>%
    tidyr::pivot_wider(names_from = percentile) %>%
    ggplot2::ggplot(ggplot2::aes(x = year, y = `50%`)) +
    #  ggplot2::scale_fill_brewer(direction = -1) +
    #ggplot2::aes_string(x = results_fits[[1]]$year, y = "`50%`")) +
    # ggplot2::geom_line(aes(y = .data$`2.5%`, color = model), lty = 2) +
    # ggplot2::geom_line(aes(y = .data$`97.5%`, color = model), lty = 2) +
    # ggplot2::geom_line(aes(y = .data$`10%`, color = model), lty = 3) +
    # ggplot2::geom_line(aes(y = .data$`90%`, color = model), lty = 3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = model), alpha = 0.1) +
    # ggplot2::geom_ribbon(aes(ymin = `10%`,  ymax = `90%`, fill = model), alpha = 0.1) +
    #   ggplot2::geom_ribbon(aes(ymin = .data$`2.5%`,  ymax = .data$`97.5%`, fill = model), alpha = 0.1) +

    ggplot2::geom_line(ggplot2::aes(color = model, lty = model), linewidth = 1.2) +
    ggplot2::labs(#fill = "Model (w 80%CI)",
      x = "year") +
    ggplot2::ylab(indicator_select) +
    #    ylab("95 & 80% interval")+
    ggplot2::theme_bw() +
    ggplot2::expand_limits(y = 0) +
    # tmp: if we don't want a model legend
  #  guides(fill = "none", lty = "none", color= "none")
    ggplot2::guides(fill = "none", lty = "none")
  if (!is.null(filtered_data)){
    # if (sum(filtered_data$held_out) > 0){
    #   filtered_data <- filtered_data %>%
    #     mutate(included = ifelse(held_out ==0, "Yes", "No"))
    # }
    filtered_data <- filtered_data %>%
      dplyr::mutate(included = ifelse(held_out ==0, "Yes", "No"))
    p <- p +
      ggnewscale::new_scale_color() +
      ggplot2::geom_errorbar(data = filtered_data,
                             ggplot2::aes(y = indicator,  ymin = indicator_low,
                                 ymax = indicator_up, color = data_series_type), alpha = 0.3) +
      ggplot2::geom_point(data = filtered_data,
                          ggplot2::aes(y = indicator,
                              x = year, color = data_series_type
                              #, shape = included
                              )) +
      ggplot2::scale_colour_manual(values = cols_sourcetypes)
  }
  return(p)
}
