#' Construct population weight matrix for one aggregate unit
#'
#' Constructs a row vector of population weights for regions in a single aggregate unit.
#'
#' @param aggregation_meta_one_unit tibble with columns agg_unit_name and region_index
#' @param population_selectedyear tibble with region_index and population, needs to include all regions
#'
#' @returns row vector with weights for each region in the aggregate unit (0 if region is not included)
#' @keywords internal
#' # right now in example_creatingaggregates.R in submeta repo
construct_popweightmatrix_for_oneaggregateunit <- function(
    aggregation_meta_one_unit, # tibble with agg_unit_name and region_index
    population_selectedyear, # tibble with region_index and population, for all regions!
    name_region_index = "region_index"
) {

  # consider starting with region_index.. ordering matters here!
  pop_wt <- population_selectedyear %>% # start with population_selectedyear because it has all regions
    left_join(aggregation_meta_one_unit) %>%
    mutate(
      # set pop to 0 if region not in aggregate
      population = ifelse(is.na(agg_unit_name), 0, population),
      weight = population/sum(population)
    ) %>%
    rename(region_index = name_region_index) %>%
    arrange(region_index) %>%
    dplyr::select(region_index, weight) %>%
    tidyr::pivot_wider(names_from = region_index, values_from = weight) %>%
    as.matrix()
  # returns row vector
  return(pop_wt)
}

