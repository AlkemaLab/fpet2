
#' Plot subnational comparison
#'
#' @param results output from fit_fpem()
#' @param indicator_select choice out of c("contraceptive_use_modern", "unmet_need_modern",
#' "demand", "demand_satisfied_modern")
#' @param year_select year to plot
#' @param marital_status_select choice out of c("all", "married", "unmarried")
#' @param ymin_select default 0
#' @param ymax_select default NA
#' @param arrange_point logical, if TRUE arrange points by value of indicator_select
#'
#' @returns a ggplot object
#' @keywords internal
plot_subnational_comparison <- function(results, indicator_select, year_select,
                                        marital_status_select,
                                        ymin_select = 0,
                                        ymax_select = NA,
                                        arrange_point = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_subnational_comparison(). Please install it.", call. = FALSE)
  }
  res <- results$estimates |>
    dplyr::filter(indicator == indicator_select, year == year_select, marital_status == marital_status_select,
                  measure == "proportion") |>

    tidyr::pivot_wider(names_from = percentile, values_from = value) |>
    dplyr::mutate(
      region_code = factor(region_code) |>
        forcats::fct_relevel("National") |>
        forcats::fct_rev())
  if (arrange_point){
    res <-  res |>
      mutate(region_code = fct_reorder(region_code, mean))
  }
  res |>
    ggplot2::ggplot(ggplot2::aes(x = region_code, y = mean,
                                 # these colors don't work as red/black but differences show up :)
                                 color = ifelse(region_code == "National", "red", "black"))) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Region", y = "Percentage",
                  title = paste0(gsub("_", " ", indicator_select), ", ", year_select)) +
    ggplot2::ylim(ymin_select, ymax_select) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_flip()
}

