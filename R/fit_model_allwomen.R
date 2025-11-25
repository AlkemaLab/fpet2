

#' Fit model to produce all women estimates
#'
#' This function is called from fit_fpem,
#' to do the model fitting for all women
#'
#' @param data_allwomen list with input data for married and data for unmarried women
#' @param runstep local_national or local_subnational
#' @param population_df pop counts for each geounit-year combi
#' @param emu_list emu info from `service_statistics_getstandata`
#' @param chains Number of chains for Stan
#' @param iter_sampling Number of iterations for sampling
#' @param iter_warmup Number of iterations for warmup
#' @param seed Random seed for Stan
#' @param add_inits Use function to sample initial values? Defaults to FALSE
#' @param add_sample TRUE
#' @param compile_model TRUE
#' @param force_recompile FALSE
#' @param refresh 10
#' @param adapt_delta 0.9
#' @param max_treedepth 14
#' @param variational FALSE
#' @param nthreads_variational 8
#'
#' @returns results list with stan_model, stan_data, samples, runstep
#' @keywords internal

fit_model_allwomen <- function(
  data_allwomen, #  this contains list with data from marital groups
  runstep = "local_national", # local_national or local_subnational
  population_df, # pop counts for each geounit-year combi
  emu_list = NULL,

  end_year,
  # settings for sampling
  chains = 4,
  iter_sampling = 200,
  iter_warmup = 150,

  add_inits = FALSE,
  add_sample = TRUE, # if FALSE, returns fit w/o samples
  # minor to do: check settings when finalizing stan model block/decision re instantiate
  compile_model  = TRUE, force_recompile = FALSE,
  seed = 1234,
  refresh = 50,
  adapt_delta = 0.9,
  max_treedepth = 14,
  # # settings for variational inference
  # not tested yet
  variational = FALSE,
  # if variational, just compile model with threading support and pass back model and stan_data
  ##variational = variational | !add_sample, # no sampling when TRUE
  nthreads_variational = 8 #40, # 8
  # max_lbfgs_iters = 1000, # default is 1000
  # num_psis_draws = 1000,



) {

  # check that geo_unit_index is ok
  if (runstep == "local_subnational"){
    checkregions <- all.equal.character(
      data_allwomen[["married"]]$geo_unit_index %>%
        arrange(c) %>%
        pull(region_code),
      data_allwomen[["unmarried"]]$geo_unit_index %>%
        arrange(c) %>%
        pull(region_code))
    if (checkregions != TRUE){
      stop("geo_unit_index not the same for married and unmarried")
    }
  }
  result <- list()
  #### stan_data
  # we use t_min from married (t_max too but that's already fixed to 2024 anyway)
  if (any(data_allwomen[["married"]]$stan_data$t_min >
          data_allwomen[["unmarried"]]$stan_data$t_min)){
    data_allwomen[["married"]]$stan_data$t_min <- rep(min(data_allwomen[["unmarried"]]$stan_data$t_min),
                                                      data_allwomen[["unmarried"]]$stan_data$n_geounit)
  }
  names(data_allwomen[["unmarried"]]$stan_data) <- paste0("unmarried_", names(data_allwomen[["unmarried"]]$stan_data))
  stan_data <- c(
    data_allwomen[["married"]]$stan_data,
    data_allwomen[["unmarried"]]$stan_data,
    is_aggr = data_allwomen[["married"]]$data$is_aggr,
    unmarried_is_aggr = data_allwomen[["unmarried"]]$data$is_aggr)

  # add population weights
  asmatrix_prop_married_t <- population_df %>%
    dplyr::group_by(year, is_in_union) %>%
    dplyr::summarise(popcount = sum(population_count)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(prop_married = sum(popcount[is_in_union == "Y"]) /
                sum(popcount[is_in_union == "Y"] + popcount[is_in_union == "N"])) %>%
    dplyr::filter(year <= end_year) %>%
    dplyr::arrange(year) %>%
    dplyr::select(-year) %>%
    as.matrix()
  if (runstep == "local_subnational"){
    stan_data$prop_married_rt <- population_df %>%
      group_by(region_code, year, is_in_union) %>%
      summarise(popcount = sum(population_count)) %>%
      group_by(region_code, year) %>%
      summarise(prop_married = sum(popcount[is_in_union == "Y"]) /
                  sum(popcount[is_in_union == "Y"] + popcount[is_in_union == "N"])) %>%
      ungroup() %>%
      arrange(year) %>%
      left_join(data_allwomen[[1]]$geo_unit_index %>%
                  dplyr::select(region_code, c)) %>%
      arrange(c) %>%
      tidyr::pivot_wider(names_from = year, values_from = prop_married) %>%
      dplyr::select(-c, -region_code)  %>%
      as.matrix()
  } else {
   #at national level, same as prop_married_rt for national run
    stan_data$prop_married_rt <- t(asmatrix_prop_married_t)
  }
  stan_data$prop_married_t <- c(asmatrix_prop_married_t)

  if (!is.null(emu_list)){
    stan_data <- c(stan_data, emu_list)
  }
  #return(stan_data)

  #### reading stan model
  if (runstep == "local_national"){
    if (!is.null(emu_list)){
      #stan_file_path <- system.file("stan/fpem_buildingblocks_complete_wtrad_allwomen_emu.stan", package = "fpet2")
      stan_file_path <- system.file("stan/fpem_allwomen_emu.stan", package = "fpet2")
    } else {
      #stan_file_path <- system.file("stan/fpem_buildingblocks_complete_wtrad_allwomen.stan", package = "fpet2")
      stan_file_path <- system.file("stan/fpem_allwomen.stan", package = "fpet2")
    }
  } else {
    if (!is.null(emu_list)){
      # stan_file_path <- system.file("stan/fpem_buildingblocks_complete_allwomenaggregates_emu.stan", package = "fpet2")
      stan_file_path <- system.file("stan/fpem_allwomen_aggregates_emu.stan", package = "fpet2")
    } else {
      #stan_file_path <- system.file("stan/fpem_buildingblocks_complete_allwomenaggregates.stan", package = "fpet2")
      stan_file_path <- system.file("stan/fpem_allwomen_aggregates.stan", package = "fpet2")
    }
  }

  ##### Load model #####
  if (compile_model & add_sample){
    stan_model <- compile_model(variational = variational,
                                nthreads_variational = nthreads_variational,
                                force_recompile = force_recompile,
                                stan_file_path = stan_file_path )
    result$stan_model <- stan_model
  }

  ##### Fit model ########
  if (add_inits){
    init_ll <- lapply(1:chains, function(id) init_fun(chain_id = id, stan_data))
  } else {
    init_ll <- NULL
  }
  #options(cmdstanr_warn_inits = FALSE)
  fit <- stan_model$sample(
    stan_data,
    save_latent_dynamics = TRUE,
    init = init_ll,
    chains = chains,
    parallel_chains  = chains,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    seed = seed,
    refresh = refresh,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    show_exceptions = FALSE
  )

  result <- c(result, list(
                stan_data = stan_data,
                samples = fit, runstep = runstep))

  return(result)

}

