

#' Fit FPET models for one marital group
#'
#' @param survey_df tibble with survey data
#' @param routine_df tibble with routine data
#' @param year column name of outcome year.
#' @param source column name of data source.
#' @param area column name of the area of each observation
#'
#' @param population_data a data frame with yearly population counts for
#'   subnational regions. It should have columns matching the names specified
#'   for \code{year} and \code{area}. This data frame is only required if the
#'   primary \code{data} set contains a mix of observations at national and
#'   subnational levels.
#'
#' @param runstep type of run (see Details).
#' @param global_fit optional object of class `"fpemplus"`, used to obtain fixed
#'   values to use for some parameters in the current fit (see Details).
#' @param iso_select ISO code to use for local runs
#'
#' @param start_year start year of estimates.
#' @param end_year end year of estimates.
#'
#' @param t_star reference year used in model.
#' @param num_knots number of spline knots.
#' @param spline_degree spline degree. Degree 2 or 3 is supported.
#' @param hierarchical_asymptote vector specifying hierarchical structure for
#'   asymptote (see Details).
#' @param hierarchical_level vector specifying hierarchical structure for the
#'   level in reference year (see Details).
#' @param hierarchical_splines vector specifying hierarchical structure for
#'   spline coefficients (see Details).
#' @param Betas_upper_bound upper bound for the splines parameters
#' @param Betas_lower_bound lower bound for the splines parameters
#' @param Ptilde_low lower bound for the asymptote Ptilde
#' @param add_dataoutliers boolean indicator of whether to include data outliers in 1b
#'
#' @param add_shocks boolean indicator of whether to include shocks (not yet included)
#' @param extra_stan_data list of additional data to pass to Stan model

#' @param get_posteriors boolean indicator of whether to return posterior samples
#'
#' @param held_out binary vector indicating which observations are held out. Set to FALSE to hold out no observations.
#' @param validation_cutoff_year year to use for out-of-sample validation, overwrites held_out (to confirm it does)
#' @param validation_run boolean indicator of whether it's a validation model run or not
#'
#' @param generate_quantities binary vector indicating whether to simulate data from the fitted model
#'
#' Setting for where to save things
#' @param output_dir output directory, defaults to NULL
#' @param runname name to use for run, if output_dir is NULL
#' @param runnumber number to add to runname
#' @param rungroup group to add to runname
#'
#' Settings for sampling
#' @param add_inits boolean indicator of whether to add initial values to the Stan model
#' @param chains number of chains to run
#' @param iter_sampling number of posterior samples to draw
#' @param iter_warmup number of warmup iterations
#' @param add_sample boolean indicator of whether to return samples
#' @param compile_model boolean indicator of whether to compile the Stan model
#' @param force_recompile boolean indicator of whether to force recompilation of the Stan model
#' @param seed random seed
#' @param refresh number of iterations between progress updates
#' @param adapt_delta target acceptance rate for the No-U-Turn Sampler
#' @param max_treedepth maximum tree depth for the No-U-Turn Sampler
#' @param variational boolean indicator of whether to use variational inference (not yet tested)
#' @param nthreads_variational number of threads to use for variational inference
#'
#' @return fpemplus object.
#'
#'
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @importFrom tibble tibble
#' @importFrom splines bs
#' @import dplyr
#' @importFrom readr read_file
#' @importFrom stringr str_replace_all
#' @import localhierarchy
#'
#' @export
#'
#'
#'
#'
fit_model <- function(
  survey_df,
  routine_df = NULL,
  population_data = NULL,

  DM1_y = "logit_contraceptive_use_modern",
  DM1_se = "se_logit_modern",
  DM2_y = "logit_unmet_overnotmodern",
  DM2_se = "se_logit_unmet_overnotmodern",

  year = "year",
  source = "data_series_type",
  area = "iso", # "iso" for national level, "region_code" for subnat

  iso_select  = NULL, # required for local national run and local_subnat run

  # type of run is defined by runstep:
  runstep, # type of run, step or localnat or localsubnat
  # step1a =  get subcluster parameters, fit w/o ar
  # step1b =  fix subcluster info and get ar and data outlier parameters and across country sigmas
  # local_national =  get local results only
  # global_subnational, local_subnational
  global_fit = NULL, # eventually, read in from data_raw if needed but NULL

  ## arguments that are relevant for first global run only (follow from global_fit in other cases)
  year_star = 2004, # for fp, 2004
          # year_star needs to be within the estimation period
  spline_degree = 2,
  num_knots = 5, # 8 for cd

  add_dataoutliers = TRUE,
  is_married = TRUE,
  #marital_status = "married", # derived from is_married

  # june 13, 2025: add back end_year as an argument
  end_year = 2030,

  # hierarchical parameters
  # used only if global_fit is NULL
  # eta so ds
  transitionmodelparam = list(
    prefix = "",
    hierarchical_terms = list(
      married = list(
        hierarchical_level  = c("intercept", "cluster", "subcluster", "iso"),
        hierarchical_splines  = c("intercept", "cluster", "subcluster", "iso"),
        hierarchical_asymptote = c("intercept", "cluster", "iso")), #"subcluster",
      # may 19
      unmarried = list(
        hierarchical_level     = c("intercept", "regional2_unmarried", "iso"),
        hierarchical_splines   = c("intercept", "clusterandsa0_unmarried", "iso"),
        hierarchical_asymptote = c("intercept",  "clusterandsa0_unmarried", "iso"))
      # # may 18
      # unmarried = list(
      #   hierarchical_level     = c("intercept", "regional2_unmarried", "iso"),
      #   hierarchical_splines   = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried", "iso"),
      #   hierarchical_asymptote = c("intercept",  "clusterandsa0_unmarried", "iso"))
    # when running just with SA1
      # unmarried = list(
      #   hierarchical_level  = c("intercept", "cluster", "subcluster", "iso"),
      #   hierarchical_splines  = c("intercept", "cluster", "subcluster", "iso"),
      #   hierarchical_asymptote = c("intercept", "cluster", "iso"))
    ),
    stan_data_settings = list(
      married = list(
        Betas_upper_bound = 0.5,
        Betas_lower_bound = 0.01,
        Ptilde_low = 0.1,
        Ptilde_scalarprior_sd = 3,
        Ptilde_scalarprior_mean = 0,
        Ptilde_prior_sd_sigma_estimate = 2,
        Omega_scalarprior_sd = 3,
        Omega_scalarprior_mean = 0,
        Omega_prior_sd_sigma_estimate = 2,
        Betas_scalarprior_sd = 3,
        Betas_scalarprior_mean = -1,
        Betas_prior_sd_sigma_estimate = 2),
      unmarried = list(
        Betas_upper_bound = 0.5,
        Betas_lower_bound = 0.01,
        Ptilde_low = 0.05,
        Ptilde_scalarprior_sd = 3,
        Ptilde_scalarprior_mean = 0,
        Ptilde_prior_sd_sigma_estimate = 2,
        Omega_scalarprior_sd = 3,
        Omega_scalarprior_mean = 0,
        Omega_prior_sd_sigma_estimate = 2,
        Betas_scalarprior_sd = 3,
        Betas_scalarprior_mean = -1,
        Betas_prior_sd_sigma_estimate = 2)
      )),

  d_transitionmodelparam = list(
    prefix = "d_",
    hierarchical_terms = list(
      married = list(
        hierarchical_level  = c("intercept", "cluster", "subcluster", "iso"),
        hierarchical_splines  = c("intercept", "cluster", "subcluster", "iso"),
        hierarchical_asymptote = c("intercept", "cluster", "iso")),
      unmarried = list(
        hierarchical_level     = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried",  "iso"),
        hierarchical_splines   = c("intercept", "clusterandsa0_unmarried",  "iso"),
        hierarchical_asymptote = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried",  "iso"))
      # may 20 v2
      # unmarried = list(
      #   hierarchical_level     = c("intercept",  "is_unmarried_sexual_activity", "regional2_unmarried",  "iso"),
      #   hierarchical_splines   = c("intercept", "clusterandsa0_unmarried",  "iso"),
      #   hierarchical_asymptote = c("intercept", "is_unmarried_sexual_activity", "clusterandsa0_unmarried",  "iso"))
      # unmarried = list(
      #   hierarchical_level     = c("intercept", "is_unmarried_sexual_activity", "clusterandsa0_unmarried",  "iso"),
      #   hierarchical_splines   = c("intercept", "clusterandsa0_unmarried",  "iso"),
      #   hierarchical_asymptote = c("intercept", "is_unmarried_sexual_activity", "clusterandsa0_unmarried",  "iso"))
      # # may 19 v1
      # unmarried = list(
      #   hierarchical_level     = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried", "iso"),
      #   hierarchical_splines   = c("intercept", "clusterandsa0_unmarried",  "iso"),
      #   hierarchical_asymptote = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried", "iso"))
      # may 18
      # unmarried = list(
      #   hierarchical_level     = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried", "iso"),
      #   hierarchical_splines   = c("intercept", "is_unmarried_sexual_activity", "hier_regional_unmarried", "iso"),
      #   hierarchical_asymptote = c("intercept", "clusterandsa0_unmarried",  "iso"))#"hier_regional_unmarried", "iso"))
      # when running just with SA1
      # unmarried = list(
      #   hierarchical_level  = c("intercept", "cluster", "subcluster", "iso"),
      #   hierarchical_splines  = c("intercept", "cluster", "subcluster", "iso"),
      #   hierarchical_asymptote = c("intercept", "cluster", "iso"))
    ),
    stan_data_settings = list(
      married = list(
        Betas_upper_bound = 0.5,
        Betas_lower_bound = 0.01,
        Ptilde_low = 0.1,
        Ptilde_scalarprior_sd = 3,
        Ptilde_scalarprior_mean = 0,
        Ptilde_prior_sd_sigma_estimate = 2,
        Omega_scalarprior_sd = 3,
        Omega_scalarprior_mean = 0,
        Omega_prior_sd_sigma_estimate = 2,
        Betas_scalarprior_sd = 3,
        Betas_scalarprior_mean = -1,
        Betas_prior_sd_sigma_estimate = 2),
      unmarried = list(
        Betas_upper_bound = 0.5,
        Betas_lower_bound = 0.01,
        Ptilde_low = 0.01,
        Ptilde_scalarprior_sd = 3,
        Ptilde_scalarprior_mean = 0,
        Ptilde_prior_sd_sigma_estimate = 2,
        Omega_scalarprior_sd = 3,
        Omega_scalarprior_mean = 0,
        Omega_prior_sd_sigma_estimate = 2,
        Betas_scalarprior_sd = 3,
        Betas_scalarprior_mean = -1,
        Betas_prior_sd_sigma_estimate = 2)
    )),

  z_transitionmodelparam = list(
    prefix = "z_",
    hierarchical_terms = list(
      married = list(
        hierarchical_level  = c("intercept", "cluster", "subcluster", "iso")),
      unmarried = list(
        hierarchical_level     = c("intercept", "clusterandsa0_unmarried", "iso"))
        ),
    stan_data_settings = list(
      married = list(
        Omega_scalarprior_sd = 3,
        Omega_scalarprior_mean = 0,
        Omega_prior_sd_sigma_estimate = 2),
      unmarried = list(
        Omega_scalarprior_sd = 3,
        Omega_scalarprior_mean = 0,
        Omega_prior_sd_sigma_estimate = 2)
    )),

  # for subnational, one level is added
  add_subnational_hierarchy = "region_code", #admin1", # this is what's added to the hierarchy for subnational
  add_aggregates  = FALSE, # to add national level aggregates, only for subnational local run with all regions included
  population_df = NULL,
  add_shocks =  FALSE, # not yet implemented
  # note: if we want to add many more arguments, can also consider adding to extra_stan_data
  extra_stan_data = list(),

  # Out-of-sample validation
  # to do: check that combi of held_out and validation_cutoff_year are still used correctly
  held_out = FALSE,
  validation_cutoff_year = NULL, # if not NULL, should be a year and is used to define/overwrite held_out set
  validation_run = FALSE,
  # Model checks
  generate_quantities = TRUE,
  # misc
  get_posteriors = TRUE,

  # outputdir
  output_dir = NULL,
  ## minor to do: this is automated, consider updating default
  #create_runname_and_outputdir = TRUE,
  runnumber = 1, # used if create_runname_and_outputdir , added to runname, increased automatically if directory exists
  rungroup = NULL,
  runname = NULL, # uswd if !create_runname_and_outputdir

  # settings for sampling
  chains = 4, # probably need more for final model
  iter_sampling = 200,
  iter_warmup = 150,
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
  nthreads_variational = 8, #40, # 8
  # max_lbfgs_iters = 1000, # default is 1000
  # num_psis_draws = 1000,

  add_inits = TRUE
  # # Stan settings
  # ...

) {


  marital_status = ifelse(is_married, "married", "unmarried")

  # to do: update naming here, here subnational refers to
  # whether or not correlations are introduced
  subnational = FALSE
  correlated_smoothing = FALSE
  correlated_smoothing_group = "iso"
  fix_subnat_corr = FALSE
  # subnational
  if (runstep %in% c("local_subnational", "global_subnational")){
    #print("TMP solution: we do a subnational run but we set subational to false to not deal with aggregates ")
    subnational = FALSE #TRUE
    print("We use subnational data.")
    # already in argument!
    #area = "admin1" #"region_code"
    #print("We do NOT add subnational correlation.")
    correlated_smoothing = FALSE #TRUE
    # if (runstep == "step3"){
    #   print("We estimate subnational correlation.")
    #   fix_subnat_corr = FALSE
    # }
    #  if (runstep %in% c("local_subnational")){
    #    print("We fix subnat correlation.")
    fix_subnat_corr = FALSE # not fixing else we look for it in global fit
    # }
  }

  data <- survey_df # for now, just to keep the same name as in fpem

  # years to produce estimates for
  # start_year needs to stay fixed for t_star and t_max so move out of arguments right now
  start_year = 1970
  if (end_year < 2030) stop("end_year less than 2030 not tested")
  # to do: tstar currently depends on keeping start_year and end_year fixed
  # if (runstep %in% c("step1a", "step1b", "step2", "global_subnational")){
  #   start_year = round(min(data$year, na.rm = TRUE)) -1
  #   end_year = round(max(data$year, na.rm = TRUE)) + 1
  # }

  indicator <- data %>% pull(indic) %>% unique()
  if(nrow(data) == 0) {
    stop("Data has no rows.")
  }

  if(start_year > end_year) {
    stop("start_year must be less than end year")
  }

  if (!runstep %in% c("step1a", "step1b", "step2", "local_national", "global_subnational", "local_subnational")){
    stop("runstep not yet implemented!")
  }

  if (! (runstep %in% c("step1a")) & is.null(global_fit)){
    stop("Need a global fit for this run.")
    # we can also read in an internal object but for now, require it to be passed
    # globalstepname <- dplyr::case_when(
    #   runstep == "step1b" ~ "1a",
    #   runstep == "local_national" ~ "1b",
    #   runstep == "step3" ~ "1b",
    #   TRUE ~ "3"
    # )
    # global_fit <- readRDS(file = paste0(
    #   here::here("data-raw/internal/"), indicator, "_summary",
    #   globalstepname,
    #   ".rds"))
  }
  if (runstep %in% c("step1a")){
    # first global fit, need to use arguments from function call
    # let's check those
    if(!(year_star %in% (start_year +1): (end_year -1)))
      stop(glue::glue("Supplied year_star {year_star} not in range of estimation years ({start_year+1} to {end_year-1})."))
    t_star <- which(seq(start_year, end_year) == year_star)
    stopifnot(is.numeric(spline_degree))
    stopifnot(is.numeric(num_knots))
    if(!(spline_degree %in% c(2, 3))) {
      stop("spline_degree must be either 2 or 3.")
    }
    if(num_knots <= 0) {
      stop("num_knots must be greater than zero.")
    }
  } else {
    # all other runs, NOT 1a
    print("We use a global fit, and take selected settings from there.")
    print("settings for the spline_degree and num_knots taken from global run")
    spline_degree <- global_fit$spline_degree
    num_knots <- global_fit$num_knots
    print("Setting for tstar taken from global run")
    t_star <- global_fit$t_star
  }
  # hierarchical levels and what's fixed also follows from run_step
  if (runstep %in% c("step1a")){
    # check settings from arguments
    if(length(transitionmodelparam$hierarchical_terms[[marital_status]]$hierarchical_asymptote) == 0) {
      stop("No hierarchical structure supplied for the asymptote. See the hierarchical_asymptote argument.")
    }
    if(length(transitionmodelparam$hierarchical_terms[[marital_status]]$hierarchical_level) == 0) {
      stop("No hierarchical structure supplied for the level in reference year. See the hierarchical_level argument.")
    }
    if(length(transitionmodelparam$hierarchical_terms[[marital_status]]$hierarchical_splines) == 0) {
      stop("No hierarchical structure supplied for the spline coefficients. See the hierarchical_splines argument.")
    }
  }
  hierarchical_terms_and_fixed <- get_hierlevels(runstep,
                                hierarchical_terms = transitionmodelparam$hierarchical_terms[[marital_status]],
                                global_fit = global_fit,
                                prefix = "",
                                add_subnational_hierarchy = add_subnational_hierarchy)
  d_hierarchical_terms_and_fixed <- get_hierlevels(runstep,
                                                 hierarchical_terms = d_transitionmodelparam$hierarchical_terms[[marital_status]],
                                                 global_fit = global_fit,
                                                 prefix = "d_",
                                                 add_subnational_hierarchy = add_subnational_hierarchy)

  smoothing <- ifelse(runstep != "step1a", TRUE, FALSE)
  # non-fixing of trad is dealt with separately
  fix_smoothing <- ifelse(runstep %in% c("step1a", "step1b"), FALSE, TRUE)
  fix_nonse <- ifelse(runstep %in% c("step1a", "step1b"), FALSE, TRUE)
  # dataoutliers is based on argument in 1a, else global fit
  add_dataoutliers <- ifelse(runstep == "step1a", TRUE, global_fit$add_dataoutliers)

  # get the spline data (needed for details on hier parameters)
  stan_spline_data <- get_stan_spline_data(num_knots = num_knots, spline_degree = spline_degree)


  ##### Data processing  #####
  # keep the names of the original data
  names_original_data <- names(data)

  # to do: just filter, warning if no data!
  # Make sure the observed data are within the estimation period
  if(sum(!(data[[year]] %in% start_year:end_year)) > 0) {
    print(glue::glue("Observations included in dataset that fall outside the estimation period ({start_year} to {end_year})."))
    print("We filter data to be inside estimation period")
    data <- data %>%
      dplyr::filter(year >= start_year, year <= end_year)
  }
  if (runstep  %in%  c("local_national", "local_subnational")){
    print(paste("Local run, we only use data for", iso_select))
    data <- data %>%
      dplyr::filter(iso %in% iso_select)
    if (dim(data)[1] == 0){
      stop("No data for this iso in the dataset")
    }
  }

  if (runstep == "step1a"){
    print("We do give nonSE to DHS (by temporarily renaming DHS into DHS0)")
    data <- data %>%
      mutate(data_series_type = ifelse(data_series_type == "DHS", "DHS0", data_series_type))
  }

  # what validation data to use?
  if (!is.null(validation_cutoff_year)){
    print("Setting held_out based on validation_cutoff_year")
    held_out <- as.numeric(data[["start_date"]] >= validation_cutoff_year)
    validation_run <- TRUE
    print(held_out)
  }
  if(length(held_out) == 1 && held_out == FALSE) {
    held_out = rep(0, nrow(data))
  }  else {
    if(length(held_out) != nrow(data)) stop(glue::glue("held_out (length {length(held_out)}) must be same size as dataset ({nrow(data)} rows)."))
    held_out = as.numeric(held_out)
  }
  data$held_out <- held_out


  ##### Setup data for Stan #####

  # Create year lookup table
  time_index <- tibble::tibble(
    year = seq(start_year, end_year),
    t = 1:length(year)
  )

  source_index <- data %>%
    dplyr::distinct("{source}" := .data[[source]]) %>%
    dplyr::mutate(source = 1:n())


  # Create district index for matching district and district index
  hierarchical_column_names <- unique(c(
    hierarchical_terms_and_fixed$hierarchical_asymptote,
    hierarchical_terms_and_fixed$hierarchical_splines,
    hierarchical_terms_and_fixed$hierarchical_level,
    d_hierarchical_terms_and_fixed$hierarchical_asymptote,
    d_hierarchical_terms_and_fixed$hierarchical_splines,
    d_hierarchical_terms_and_fixed$hierarchical_level
  )) %>%
    setdiff("intercept")

  # Make sure there are no NAs in any of the columns
  for(column in hierarchical_column_names) {
    if(column == "intercept") next
    #if (column == area) {
    #  check_nas_or_pops(data, column, year, population_data)
    #} else {
    localhierarchy::check_nas(data, column)
    #}
  }

  # # update 5/23: if national data was included, this is excluded when creating the geo_unit_index
  # set regio_code to NA if national
  #print(data[[area]])
  data[[area]] <- ifelse(data[[area]] == "National", NA, data[[area]])
  #print(data[[area]])
  geo_unit_index <- data[!is.na(data[[area]]), ] %>%
    dplyr::distinct(!!! syms(hierarchical_column_names)) %>%
    dplyr::mutate(c = 1:n())



  year_by <- c()
  year_by[year] = "year"
  data <- data %>%
    dplyr::left_join(time_index, by = year_by) %>%
    dplyr::left_join(geo_unit_index, by = hierarchical_column_names) %>%
    dplyr::left_join(source_index, by = source)

  # to do: make sure this is not a requirement
  # for fp, data needs to be sorted by country
  data <- data %>% arrange(c)
  # now save original data
  original_data <- data %>%
    dplyr::select(all_of(names_original_data))


  ## all data checks and imputation related to NAs
  # Make sure there are no NAs in supplied columns
  check_nas(data, year)
  check_nas(data, source)

  # maybe less relevant for 1 indicator...
  # create obs_isNA
  DM1_obs_isna <- is.na(data[[DM1_y]])
  DM2_obs_isna <- is.na(data[[DM2_y]])
  check_nas(data[!DM1_obs_isna, ], DM1_se)
  check_nas(data[!DM2_obs_isna, ], DM2_se)
  # impute arbitrary numbers to avoid NA issues?
  # note: now we can't use the DM1_y etc anymore below!!!
  # to update to skip over in a next version
  data[[DM1_y]][DM1_obs_isna] <- 100
  data[[DM1_se]][DM1_obs_isna] <- 100
  data[[DM2_y]][DM2_obs_isna] <- 100
  data[[DM2_se]][DM2_obs_isna] <- 100

  # Make sure SEs are all positive and non-zero
  if(sum(data[[DM1_se]] <= 0) > 0) {
    stop(glue::glue("All standard errors must be greater than zero ({sum(data[[DM1_se]] <= 0)} observations supplied with zero or negative SEs)."))
  }
  if(sum(data[[DM2_se]] <= 0) > 0) {
    stop(glue::glue("All standard errors must be greater than zero ({sum(data[[DM2_se]] <= 0)} observations supplied with zero or negative SEs)."))
  }

  # In "local" fits, ensure that for all hierarchy levels where any quantity is
  # fixed, all geographic units in data used for the local fit also were present
  #  in the global fit

  # checks
  # june 24, remove check
  # if (!is.null(global_fit)) {
  #   fixed_hierarchy_levels <- unique(c(hierarchical_terms_and_fixed$hierarchical_asymptote_sigmas_fixed,
  #                                      hierarchical_terms_and_fixed$hierarchical_asymptote_terms_fixed,
  #                                      hierarchical_terms_and_fixed$hierarchical_splines_sigmas_fixed,
  #                                      hierarchical_terms_and_fixed$hierarchical_splines_terms_fixed,
  #                                      hierarchical_terms_and_fixed$hierarchical_level_sigmas_fixed,
  #                                      hierarchical_terms_and_fixed$hierarchical_level_terms_fixed,
  #                                      d_hierarchical_terms_and_fixed$hierarchical_asymptote_sigmas_fixed,
  #                                        d_hierarchical_terms_and_fixed$hierarchical_asymptote_terms_fixed,
  #                                        d_hierarchical_terms_and_fixed$hierarchical_splines_sigmas_fixed,
  #                                        d_hierarchical_terms_and_fixed$hierarchical_splines_terms_fixed,
  #                                        d_hierarchical_terms_and_fixed$hierarchical_level_sigmas_fixed,
  #                                        d_hierarchical_terms_and_fixed$hierarchical_level_terms_fixed))
  #   fixed_hierarchy_levels <- fixed_hierarchy_levels[fixed_hierarchy_levels != "intercept"]
  #
  #   fixed_geo_unit_index_local <- geo_unit_index[fixed_hierarchy_levels] |>
  #     dplyr::distinct()
  #   fixed_geo_unit_index_global <- global_fit$geo_unit_index[fixed_hierarchy_levels] |>
  #     dplyr::distinct()
  #   missing_geos <- fixed_geo_unit_index_local |>
  #     dplyr::anti_join(fixed_geo_unit_index_global,
  #                      by = fixed_hierarchy_levels)
  #   if (nrow(missing_geos) > 0) {
  #     # to do
  #     print("All geographic units that appear in the data for the current fit at hierarchical levels for which any parameter is fixed must have also been included in the data used for the `global_fit`.")
  #   }
  # }

  data <- data %>%
    # add index for DHS
    # minor to do: remove dependencies on data_series_type?
    mutate(isDHS = ifelse(data_series_type == "DHS", 1, 0))
  # to do: check if run w/o dataoutliers
  if (!add_dataoutliers){
    data <- data %>% mutate(nooutlier = 1)
  }

  # add obs period per area
  # exception is if there is national data included
  # then t_min is repeated as the min of all data, across geo units
  if (any(is.na(data[[area]]))){
    t_min <- rep(min(data$t), nrow(geo_unit_index))
  } else {
    t_min <- data %>%
      group_by(c) %>%
      summarise(t_min = min(t)) %>%
      pull(t_min)
    t_min <- ifelse(t_min > (t_star-1), t_star- 1, t_min)
  }
  # t_max <- data %>%
  #   group_by(c) %>%
  #   summarise(t_max = max(t)) %>%
  #   pull(t_max)
  # t_max <- ifelse(t_max < (t_star+ 1), t_star+ 1, t_max)
  # update for webtool runs: setting t_max to 2024 to avoid issues when including EMUs
  # this includes same dependence on start and end year as does t star
  t_max <- rep(which(seq(start_year, end_year) == 2024), nrow(geo_unit_index))

  ### PMA data processing fun
  data <- data %>%
    mutate(ispma = ifelse(data_series_type == "PMA", 1, 0),
           exact_time = (start_date + end_date)/2)
  # data are sorted by country earlier on (arrange(c))
  # for data model with correlation in pma data
  # to do: aggregate stuff
  pma_country_indices <- unique(data %>%
                                  dplyr::filter(ispma == 1 &
                                           # currently not using aggregate PMA data
                                           # is_agg_obs_tmp == 0 &
                                           # this is to create training data, remove held out data
                                           held_out ==0
                                  ) %>%
                                  pull(c))
  indices_pma <- which(data$ispma == 1 #& data$is_agg_obs_tmp == 0
                       & data$held_out == 0)
  nstart_pma_c <- npma_c <- rep(0, nrow(geo_unit_index)) #easy to make this have country indices C) #lenght(pma_country_indices))
  index <- 0
  for (index_country in pma_country_indices){
    nstart_pma_c[index_country] <- index +1
    indices_inc = which(data$c == index_country & data$ispma == 1
                        #& data$is_agg_obs_tmp == 0 #) #indices of pma data for this country
                        & data$held_out == 0)
    npma_c[index_country] = length(indices_inc)
    index <- index + npma_c[index_country]
  }

  # # check that pma doesn't have missing trad and unmet_modern entries for ALL steps
  # tmp <- data %>%
  #   dplyr::filter(data_series_type == "PMA") %>%
  #   dplyr::filter(!is.na(contraceptive_use_traditional) & is.na(unmet))
  # if (dim(tmp)[1] != 0){
  #   stop("PMA data on trad is given but PMA data on unmet is missing, we don't allow for that.")
  # }
  # if (dim(nat_data %>% dplyr::filter(data_series_type == "PMA"))[1] > 0){
  #   print("We do not allow for using national level PMA data in a subnat run. These are excluded")
  #   nat_data <- nat_data %>%
  #     dplyr::filter(data_series_type != "PMA")


  # put together stan data
  stan_data <- list(
    n_geounit = nrow(geo_unit_index),
    T = nrow(time_index),
    N = nrow(data),
    S = nrow(source_index),
    held_out = held_out,
    isdhs = data$isDHS,
    time = array(data$t),
    t_min = array(t_min),
    t_max = array(t_max),
    # for geo_unit:
    # in case of NA values in data$c, replace with the dummy value 0
    # this is only used in multiscale fitting with mixed national and subnational data
    geo_unit = array(ifelse(is.na(data$c), 0L, data$c)),
    source = array(data$source),

    DM1_s = array(data[[DM1_se]]),
    DM1_y = array(data[[DM1_y]]),
    DM1_obs_isna = DM1_obs_isna,
    # the same right now
    DM1_nooutlier = data$nooutlier,
    DM1_any_bias = data$any_bias,

    DM2_s = array(data[[DM2_se]]),
    DM2_y = array(data[[DM2_y]]),
    DM2_obs_isna = DM2_obs_isna,
    # the same right now
    DM2_nooutlier = data$nooutlier,
    DM2_any_bias = data$any_bias,

    # pma stuff
    ispma = data$ispma,
    exact_time = data$exact_time,
    pma_country_indices = pma_country_indices,
    N_pma_country_indices = length(pma_country_indices),
    npma_c  = npma_c,
    nstart_pma_c = nstart_pma_c,
    indices_pma  = indices_pma,
    N_pma_obs_indices = length(indices_pma),

    t_star = t_star,
   # t_last = t_last,
    add_dataoutliers = add_dataoutliers,
    generate_quantities = generate_quantities,
    validation_run = validation_run,
    verysmallnumber = 0.00001 # lower bound for sds
    # used earlier in generated quanities to get splines (fits)
    #eta_grid = seq(0.01, 0.99, 0.01),
    #N_eta = length(seq(0.01, 0.99, 0.01)),
  )
  #return(stan_data)


  ##### Set up hierarchical structures ######
  hier_data_transition <- get_hier_data_transition(geo_unit_index,
                                                   prefix = transitionmodelparam$prefix,
                                                   stan_data_settings = transitionmodelparam$stan_data_settings[[marital_status]],
                                                   stan_spline_data, # needed for dimension of splines coeff vector
                                                   # info on levels and what's fixed from other function
                                                   hierarchical_terms_and_fixed = hierarchical_terms_and_fixed,
                                                   global_fit = global_fit)
  d_hier_data_transition <- get_hier_data_transition(geo_unit_index,
                                                    prefix = d_transitionmodelparam$prefix,
                                                   stan_data_settings = d_transitionmodelparam$stan_data_settings[[marital_status]],
                                                   stan_spline_data, # needed for dimension of splines coeff vector
                                                   # info on levels and what's fixed from other function
                                                   hierarchical_terms_and_fixed = d_hierarchical_terms_and_fixed,
                                                   global_fit = global_fit)


  ##### Set up handling of smoothing #####
  smoothing_data <- list(
    smoothing = as.integer(smoothing),
    fix_smoothing = as.integer(fix_smoothing))
  smoothing_data <- c(smoothing_data,
                      get_smoothing_standata(fix_smoothing,
                                                         global_fit,
                                                         prefix = ""),
                      get_smoothing_standata(fix_smoothing,
                                             global_fit,
                                             prefix = "d_"))

        #### correlated stuff
        smoothing_data <- c(smoothing_data,
                            list(
          "correlated_smoothing" = as.integer(correlated_smoothing),
          "fix_subnat_corr" = as.integer(fix_subnat_corr),
          "rho_correlationeps_fixed" = numeric(0),
          "d_rho_correlationeps_fixed" = numeric(0)))
        if (fix_subnat_corr){
          if (is.null(global_fit)) {
            stop("fix_subnat_corr was set to TRUE, but a global_fit was not provided.")
          }
          smoothing_data$rho_correlationeps_fixed <-
            global_fit$post_summ %>%
            dplyr::filter(variable == "rho_correlationeps[1]") %>%
            pull(postmean)
          smoothing_data$d_rho_correlationeps_fixed <-
            global_fit$post_summ %>%
            dplyr::filter(variable == "d_rho_correlationeps[1]") %>%
            pull(postmean)
        }
        # compute info for how to do correlated smoothing
        if (correlated_smoothing) {
          if (!(correlated_smoothing_group %in% hierarchical_column_names)) {
            stop("`hierarchical_column_names` must be one of the variables used for hierarchical levels.")
          }

          # subnat corr to do: This code assumes that in the geo_unit_index,
          # lowest-level geo areas that are within the same `correlated_smoothing_group`
          # (e.g. within the same country) are assigned consecutive indices c.
          # This depends on the order of data being passed in, so maybe we should
          # also arrange the geo_unit_index by `hierarchical_column_names` at the time
          # it is created, before we use it?
          unique_group_vals <- unique(geo_unit_index[[correlated_smoothing_group]])

          smoothing_data$n_cor_smoothing_blocks <- length(unique_group_vals)
          smoothing_data$cor_smoothing_block_sizes <- array(purrr::map_dbl(
            unique_group_vals,
            function(v) {
              sum(geo_unit_index[[correlated_smoothing_group]] == v)
            }
          ))
          smoothing_data$max_cor_smoothing_block_size <- max(smoothing_data$cor_smoothing_block_sizes)
        } else {
          smoothing_data$n_cor_smoothing_blocks <- 0L
          smoothing_data$cor_smoothing_block_sizes <- integer(0)
          smoothing_data$max_cor_smoothing_block_size <- 0L
        }


  ##### Set up handing of data model hyperparameters ######
  nonse_data <- list(
      "fix_nonse" = as.integer(fix_nonse),
      DM1_fix_nonse = as.integer(fix_nonse),
      DM2_fix_nonse = as.integer(fix_nonse))

  nonse_data <- c(nonse_data,
          get_nonse_standata(fix_nonse, global_fit,
                                 source_index,
                                 prefix = "DM1_"),
          get_nonse_standata(fix_nonse, global_fit,
                             source_index,
                             prefix = "DM2_"))

  ##### Extra set up for subnational runs ######
  if (!subnational){
    #print("NOT using or producing aggregates")
    # extra stan data not used but needs to be provided
    n_agg_units <- 1
    extra_stan_data <- c(extra_stan_data,
                         list(
                           subnational = subnational,
                           is_agg_obs = rep(0, dim(data)[1]),
                           n_agg_units = n_agg_units,
                           agg_unit = rep(0, dim(data)[1]),
                           geo_unit_pop_wt = array(0, c(n_agg_units, nrow(geo_unit_index),
                                                        nrow(time_index)))
                         ))
  } else {
    # get population-based weights for each year
    # currently assumes the fit is to one aggregated (national) unit
    # this should replace the check that's done in check_nas_or_pop
    # geo_unit_subindex <- geo_unit_index[c(area, "c")]
    # required_pop_rows <- tidyr::expand_grid(
    #   year = time_index$year,
    #   area = geo_unit_subindex[[area]]) %>%
    #   dplyr::left_join(time_index, by="year")
    # colnames(required_pop_rows) <- c(year, area, "t")
    # required_pop_rows <- required_pop_rows %>%
    #   dplyr::left_join(geo_unit_subindex, by=area) %>%
    #   dplyr::left_join(population_data, by = c(year, area))
    # missing_pop_rows <- required_pop_rows |>
    #   dplyr::filter(is.na(population))
    # if (nrow(missing_pop_rows) > 0) {
    #   stop(glue::glue("If there are NAs in column {area}, `population_data` must include population data for all areas and years."))
    # }
    # geo_unit_pop_wt <- required_pop_rows[c("t", "c", "population")] %>%
    #   tidyr::pivot_wider(names_from = "t", values_from = "population") %>%
    #   dplyr::select(-c) %>%
    #   as.matrix()
    # geo_unit_pop_wt <- sweep(geo_unit_pop_wt, 2, apply(geo_unit_pop_wt, 2, sum), `/`)


    agg_data <- construct_agg_data(geo_unit_index=geo_unit_index,
                                   time_index=time_index,
                                   aggregation_meta=population_data)
    agg_geo_unit_index <- agg_data$agg_geo_unit_index
    geo_unit_pop_wt <- agg_data$geo_unit_pop_wt


    extra_stan_data <- c(extra_stan_data,


                         list(


                           subnational = subnational,
                           # indicator of whether observations are for aggregates of multiple areas,
                           # and the population weightings for those areas
                           n_agg_units = nrow(agg_geo_unit_index),
                           is_agg_obs = as.integer(data$is_agg_obs),
                           agg_unit = ifelse(
                             data$is_agg_obs,
                             data["agg_unit_name"] |>
                               dplyr::left_join(agg_geo_unit_index, by = "agg_unit_name") |>
                               dplyr::pull(i),
                             0L),
                           geo_unit_pop_wt = geo_unit_pop_wt))
  } # end subnational


  #### adding tradtional
  if (!runstep %in% c("step1a", "step1b")){
    add_trad  = TRUE
    #print("adding traditional use")
    # define the DM3_y and DM3_se columns
    # these include the exception when unmet is missing but trad is given
    # hack: input logit_trad_overnotmodern when logit_unmet_overnotmodern is na, se too
    # NEED TO USE THE DM2 INDICATORS, THE VARIABLES THEMSELVES HAVE BEEN OVERWRITTEN!
    data <- data %>%
      mutate(logit_trad_overunmet = ifelse(DM2_obs_isna, logit_trad_overnotmodern, logit_trad_overunmet),
             se_logit_trad_overunmet = ifelse(DM2_obs_isna, se_logit_trad_overnotmodern, se_logit_trad_overunmet))
    DM3_y <- "logit_trad_overunmet"
    #print(data[[DM3_y]])
    DM3_se <- "se_logit_trad_overunmet"
    DM3_obs_isna <- is.na(data[[DM3_y]])
    check_nas(data[!DM3_obs_isna, ], DM3_se)
    # impute arbitrary numbers to avoid NA issues?
    data[[DM3_y]][DM3_obs_isna] <- 100
    data[[DM3_se]][DM3_obs_isna] <- 100
    # Make sure SEs are all positive and non-zero
    if(sum(data[[DM3_se]] <= 0) > 0) {
      stop(glue::glue("All standard errors must be greater than zero ({sum(data[[DM3_se]] <= 0)} observations supplied with zero or negative SEs)."))
    }
    z_fix_smoothing <- ifelse(runstep == "step2", FALSE, TRUE)
    DM3_fix_nonse <- ifelse(runstep == "step2", FALSE, TRUE)

    stan_data <- c(stan_data,
                   list(
      DM3_s = array(data[[DM3_se]]),
      DM3_y = array(data[[DM3_y]]),
      DM3_obs_isna = DM3_obs_isna,
      DM3_nooutlier = data$nooutlier,
      DM3_any_bias = data$any_bias,
      z_fix_smoothing = z_fix_smoothing,
      DM3_fix_nonse = DM3_fix_nonse))

    # add info on parameters: z_Omega,
    z_hierarchical_terms_and_fixed <- get_hierlevels_omegatrad(runstep,
                                                     hierarchical_terms = z_transitionmodelparam$hierarchical_terms[[marital_status]],
                                                     global_fit = global_fit,
                                                     prefix = z_transitionmodelparam$prefix,
                                                     add_subnational_hierarchy = add_subnational_hierarchy)

    z_hier_data_transition <- get_hier_data_transition_justomega(geo_unit_index,
                                                       prefix = z_transitionmodelparam$prefix,
                                                       stan_data_settings = z_transitionmodelparam$stan_data_settings[[marital_status]],
                                                       # info on levels and what's fixed from other function
                                                       hierarchical_terms_and_fixed = z_hierarchical_terms_and_fixed,
                                                       global_fit = global_fit,
                                                       runstep = runstep)
    #print(z_hier_data_transition$hier_stan_data)
    # z_Rho and z_Tau for eps,
    smoothing_data <- c(smoothing_data,
                        get_smoothing_standata(fix_smoothing = z_fix_smoothing,
                                               global_fit,
                                               prefix = "z_"))
    # add DM3 model parameters
    nonse_data <- c(nonse_data,
      get_nonse_standata(fix_nonse = DM3_fix_nonse ,
                          global_fit,
                         source_index,
                         prefix = "DM3_"))


    # define model file in later block
    # end adding trad for step 2
  } else {
    add_trad = FALSE
    z_hierarchical_terms_and_fixed <- NULL
    z_hier_data_transition <- NULL
    DM3_y = NULL
    DM3_se = NULL
  }


  ### adding info for aggregates
  if (add_aggregates & runstep == "local_subnational"){
    # to do: check that population data is provided for that country only
    # check run
    # create pop weights

    #
    pop_df <- population_df %>%
      #dplyr::filter(#division_numeric_code == division_numeric_code_select,
      #       is_in_union == ifelse(is_married, "Y", "N")) %>%
      left_join(geo_unit_index) %>%
      rename(population = population_count) %>%
      left_join(time_index)

    aggregation_meta_one_unit <-
      tibble(
        agg_unit_name = rep("National", dim(geo_unit_index)[1]),
        c = geo_unit_index$c)
    popweights_tr <- matrix(NA, nrow = dim(time_index)[1],
                            ncol = dim(geo_unit_index)[1])
    for (t_select in 1:dim(time_index)[1]){
      popweights_tr[t_select,] <- construct_popweightmatrix_for_oneaggregateunit(
        aggregation_meta_one_unit = aggregation_meta_one_unit,
        population_selectedyear = pop_df %>%
          dplyr::filter(t == t_select),
        name_region_index = "c"
      )
    }
    #popweights_tr
    stan_data$geo_unit_natpop_weight_tr <- popweights_tr

  } # end adding aggregate info

  ##### Create list with combined inputs/outputs ####

  # shocks: not used yet
  extra_stan_data_shocks = list(
    add_shock = add_shocks)
  extra_stan_data <- c(extra_stan_data,
                       extra_stan_data_shocks)



  # to pass to stan
  stan_data <- c(stan_spline_data,
                 extra_stan_data,
                 stan_data,
                  hier_data_transition$hier_stan_data,
                 d_hier_data_transition$hier_stan_data,
                 z_hier_data_transition$hier_stan_data,
                 smoothing_data,
                 nonse_data)
  #return(stan_data)
  # pass back to user
  result <- list(
    runstep = runstep,

    record_id_fixed_used = data$record_id_fixed, # used to define outliers from 1a, to use in 1b

    # relabeled original data to avoid confusion
    original_data = original_data,

    # data now corresponds to input to stan_data
    data = data,
    stan_data = stan_data,
    time_index = time_index,
    geo_unit_index = geo_unit_index,
    source_index = source_index,
    year = year,
    source = source,
    area = area,

    DM1_y = DM1_y,
    DM1_se = DM1_se,
    DM2_y = DM2_y,
    DM2_se = DM2_se,
    DM3_y = DM3_y,
    DM3_se = DM3_se,

    # used in 1a only
    t_star  = t_star,
    num_knots = num_knots,
    spline_degree = spline_degree,
    smoothing = smoothing,
    add_dataoutliers = add_dataoutliers,

    held_out = held_out,
    # tau_prior = tau_prior,
    # rho_prior = rho_prior,
    # fix_subnat_corr = fix_subnat_corr,
    correlated_smoothing = correlated_smoothing,

    hierarchical_terms_and_fixed = hierarchical_terms_and_fixed,
    d_hierarchical_terms_and_fixed = d_hierarchical_terms_and_fixed,
    z_hierarchical_terms_and_fixed = z_hierarchical_terms_and_fixed
  )

  # if we want a list with hier_data
  result <- c(result, hier_data = list(hier_data_transition$hier_data, d_hier_data_transition$hier_data,
                                       z_hier_data_transition$hier_data))
  # also need this list as it's read in directly in global fit
  result <- c(result, hier_stan_data_settings = list(hier_data_transition$hier_stan_data_settings))
  result <- c(result, d_hier_stan_data_settings = list(d_hier_data_transition$hier_stan_data_settings))
  result <- c(result, z_hier_stan_data_settings = list(z_hier_data_transition$hier_stan_data_settings))
  # if we want directly in result, use append
  result <- append(result, hier_data_transition$hier_data)
  result <- append(result, d_hier_data_transition$hier_data)
  result <- append(result, z_hier_data_transition$hier_data)
  result$is_married <- is_married

  if (subnational){
    result <- c(result,
                list(agg_geo_unit_index = agg_geo_unit_index,
                     geo_unit_pop_wt = geo_unit_pop_wt))
  }



  ######
  if (!add_sample){
    return(result)
  }

  ###### reading stan model
  # here options w/o EMUs and w/o aggregates can be used (as EMUs or aggr are added only in the all women model)
  if (!add_aggregates){
    if (!add_trad){
      stan_file_path <- system.file("stan/fpem_notallwomen_notrad.stan", package = "fpet2")
    } else {
      stan_file_path <- system.file("stan/fpem_notallwomen_wtrad.stan", package = "fpet2")
    }
  } else {
    stop("We don't have a stan model for one marital group with aggregates. Did you mean to fit an all-women model with aggregates?")
  }

  ##### Load model #####
  if (compile_model & add_sample){
    stan_model <- compile_model(variational = variational,
                                nthreads_variational = nthreads_variational,
                                force_recompile = force_recompile,
                                stan_file_path = stan_file_path )
    result$stan_model <- stan_model
  }

  ##### Create an output directory for the model ####
  if (is.null(output_dir)){
    if(is.null(runname)){
      run_type <- if(validation_run == TRUE) "val" else "run"
      # set up directory to store the run
      runname <- paste0(indicator, "_", runstep, "_", run_type, "_", runnumber,
                        ifelse(variational, "_variational", ""))
      output_dir <- get_relative_output_dir(runname)
      while(dir.exists(output_dir) & runnumber < 100) {
        print("output directory already exists, increasing runnumber by 1")
        runnumber <- runnumber + 1
        runname <- paste0(indicator, "_", runstep, "_", run_type, "_", runnumber,
                          ifelse(variational, "_variational", ""))
        output_dir <- get_relative_output_dir(runname)
      }
      if (runnumber == 100){
        stop("runnumber is 100, have you really done this run 100 times already?")
      }
    } else {
      if (is.null(rungroup)) {
        output_dir <- get_relative_output_dir(runname)
      } else {
        stop("not implemented")
        #output_dir <- get_relative_output_dir(rungroup, runname)
      }
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  print(paste("output directory is", output_dir))
  result$output_dir <- output_dir

  ##### Fit model ########
  if (!variational){
    if (add_inits){
      init_ll <- lapply(1:chains, function(id) init_fun(chain_id = id, stan_data))
    } else {
      init_ll <- NULL
    }
    # can try this too
    #init = function(chain_id) init_fun(chain_id = chain_id, fit$stan_data),

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

    result <- c(result,
                samples = fit)
  } else {
    # variational
    # no longer tested
    fit_pf <- result$stan_model$pathfinder(data = result$stan_data,
                                           seed = seed,
                                        init = function(chain_id) init_fun(chain_id = chain_id,
                                                                           result$stan_data),
                                        #init = 0,
                                        # num_psis_draws = num_psis_draws,
                                        #num_paths= nthreads_variational,
                                        #num_threads = nthreads_variational,
                                        #max_lbfgs_iters = max_lbfgs_iters, # default is 1000
                                        #single_path_draws=50,
                                        output_dir = output_dir
                                        #,
                                        #history_size = 50
    )
    result <- c(result,
                samples = fit_pf)
  }



  saveRDS(result, file.path(output_dir, paste0(indicator, "_fit_nosumm.rds")))

  if (runstep %in% c("step1a", "step1b", "global_subnational", "step2")){
    result$post_summ <- get_posterior_summaries_andfindpar(result)
    saveRDS(result, file.path(output_dir, paste0(indicator, "_fit_wpostsumm.rds")))
  } #else{
  #  saveRDS(result, file.path(output_dir, paste0(indicator, "_fit_nopost.rds")))
  #}
  if (get_posteriors){
    #result$data <- add_uncertainty_in_obs(result)
    # errors for 1a/1b becasue there are no dm3 obs
    cat("Extracting posteriors...\n")
    # result$posteriors <- process_fit(result, parallel_chains = ifelse(is.null(chains), 1, chains),
    #                                  save_eps = FALSE,
    #                                  save_nontemporal  = FALSE)
    # not sure we still want this class
    # attr(result, "class") <- "fpemplus"
    result$estimates <- get_estimates_globalruns(result$samples, result$geo_unit_index,
                                                 result$time_index,
                                                 marital_status = ifelse(result$is_married, "married", "unmarried"),
                                                 return_samples = TRUE)
  }
  # if (runstep %in% c("step1a", "step1b", "global_subnational", "step2")){
  #   stepname <- dplyr::case_when(
  #     runstep == "step1a" ~ "1a",
  #     runstep == "step1b" ~ "1b",
  #     TRUE ~ runstep
  #   )
  #   summary <- result
  #   summary$samples <- NULL
  #   saveRDS(summary, file.path(output_dir, paste0(indicator, "_summary", stepname, ".rds")))
  # }

  return(result)

}

