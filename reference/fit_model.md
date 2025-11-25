# Fit FPET models for one marital group

Fit FPET models for one marital group

## Usage

``` r
fit_model(
  survey_df,
  routine_df = NULL,
  population_data = NULL,
  DM1_y = "logit_contraceptive_use_modern",
  DM1_se = "se_logit_modern",
  DM2_y = "logit_unmet_overnotmodern",
  DM2_se = "se_logit_unmet_overnotmodern",
  year = "year",
  source = "data_series_type",
  area = "iso",
  iso_select = NULL,
  runstep,
  global_fit = NULL,
  year_star = 2004,
  spline_degree = 2,
  num_knots = 5,
  add_dataoutliers = TRUE,
  is_married = TRUE,
  end_year = 2030,
  transitionmodelparam = list(prefix = "", hierarchical_terms = list(married =
    list(hierarchical_level = c("intercept", "cluster", "subcluster", "iso"),
    hierarchical_splines = c("intercept", "cluster", "subcluster", "iso"),
    hierarchical_asymptote = c("intercept", "cluster", "iso")), unmarried =
    list(hierarchical_level = c("intercept", "regional2_unmarried", "iso"),
    hierarchical_splines = c("intercept", "clusterandsa0_unmarried", "iso"),
    hierarchical_asymptote = c("intercept", "clusterandsa0_unmarried", "iso"))),
    stan_data_settings = list(married = list(Betas_upper_bound = 0.5, 
    
    Betas_lower_bound = 0.01, Ptilde_low = 0.1, Ptilde_scalarprior_sd = 3,
    Ptilde_scalarprior_mean = 0, Ptilde_prior_sd_sigma_estimate = 2, Omega_scalarprior_sd
    = 3, Omega_scalarprior_mean = 0, Omega_prior_sd_sigma_estimate = 2,
    Betas_scalarprior_sd = 3, Betas_scalarprior_mean = -1, Betas_prior_sd_sigma_estimate
    = 2), unmarried = list(Betas_upper_bound = 0.5, Betas_lower_bound = 0.01, Ptilde_low
    = 0.05, Ptilde_scalarprior_sd = 3, Ptilde_scalarprior_mean = 0,
    Ptilde_prior_sd_sigma_estimate = 2, Omega_scalarprior_sd = 3, 
    
    Omega_scalarprior_mean = 0, Omega_prior_sd_sigma_estimate = 2, Betas_scalarprior_sd =
    3, Betas_scalarprior_mean = -1, Betas_prior_sd_sigma_estimate = 2))),
  d_transitionmodelparam = list(prefix = "d_", hierarchical_terms = list(married =
    list(hierarchical_level = c("intercept", "cluster", "subcluster", "iso"),
    hierarchical_splines = c("intercept", "cluster", "subcluster", "iso"),
    hierarchical_asymptote = c("intercept", "cluster", "iso")), unmarried =
    list(hierarchical_level = c("intercept", "is_unmarried_sexual_activity",
    "hier_regional_unmarried", "iso"), hierarchical_splines = c("intercept",
    "clusterandsa0_unmarried", "iso"), hierarchical_asymptote = c("intercept",
    "is_unmarried_sexual_activity", 
     "hier_regional_unmarried", "iso"))),
    stan_data_settings = list(married = list(Betas_upper_bound = 0.5, Betas_lower_bound =
    0.01, Ptilde_low = 0.1, Ptilde_scalarprior_sd = 3, Ptilde_scalarprior_mean = 0,
    Ptilde_prior_sd_sigma_estimate = 2, Omega_scalarprior_sd = 3, Omega_scalarprior_mean
    = 0, Omega_prior_sd_sigma_estimate = 2, Betas_scalarprior_sd = 3,
    Betas_scalarprior_mean = -1, Betas_prior_sd_sigma_estimate = 2), unmarried =
    list(Betas_upper_bound = 0.5, Betas_lower_bound = 0.01, Ptilde_low = 0.01, 
    
    Ptilde_scalarprior_sd = 3, Ptilde_scalarprior_mean = 0,
    Ptilde_prior_sd_sigma_estimate = 2, Omega_scalarprior_sd = 3, Omega_scalarprior_mean
    = 0, Omega_prior_sd_sigma_estimate = 2, Betas_scalarprior_sd = 3,
    Betas_scalarprior_mean = -1, Betas_prior_sd_sigma_estimate = 2))),
  z_transitionmodelparam = list(prefix = "z_", hierarchical_terms = list(married =
    list(hierarchical_level = c("intercept", "cluster", "subcluster", "iso")), unmarried
    = list(hierarchical_level = c("intercept", "clusterandsa0_unmarried", "iso"))),
    stan_data_settings = list(married = list(Omega_scalarprior_sd = 3,
    Omega_scalarprior_mean = 0, Omega_prior_sd_sigma_estimate = 2), unmarried =
    list(Omega_scalarprior_sd = 3, Omega_scalarprior_mean = 0,
    Omega_prior_sd_sigma_estimate = 2))),
  add_subnational_hierarchy = "region_code",
  add_aggregates = FALSE,
  population_df = NULL,
  add_shocks = FALSE,
  extra_stan_data = list(),
  held_out = FALSE,
  validation_cutoff_year = NULL,
  validation_run = FALSE,
  generate_quantities = TRUE,
  get_posteriors = TRUE,
  output_dir = NULL,
  runnumber = 1,
  rungroup = NULL,
  runname = NULL,
  chains = 4,
  iter_sampling = 200,
  iter_warmup = 150,
  add_sample = TRUE,
  compile_model = TRUE,
  force_recompile = FALSE,
  seed = 1234,
  refresh = 50,
  adapt_delta = 0.9,
  max_treedepth = 14,
  variational = FALSE,
  nthreads_variational = 8,
  add_inits = TRUE
)
```

## Arguments

- survey_df:

  tibble with survey data

- routine_df:

  tibble with routine data

- population_data:

  a data frame with yearly population counts for subnational regions. It
  should have columns matching the names specified for `year` and
  `area`. This data frame is only required if the primary `data` set
  contains a mix of observations at national and subnational levels.

- year:

  column name of outcome year.

- source:

  column name of data source.

- area:

  column name of the area of each observation

- iso_select:

  ISO code to use for local runs

- runstep:

  type of run (see Details).

- global_fit:

  optional object of class \`"fpemplus"\`, used to obtain fixed values
  to use for some parameters in the current fit (see Details).

- spline_degree:

  spline degree. Degree 2 or 3 is supported.

- num_knots:

  number of spline knots.

- add_dataoutliers:

  boolean indicator of whether to include data outliers in 1b

- end_year:

  end year of estimates.

- add_shocks:

  boolean indicator of whether to include shocks (not yet included)

- extra_stan_data:

  list of additional data to pass to Stan model

- held_out:

  binary vector indicating which observations are held out. Set to FALSE
  to hold out no observations.

- validation_cutoff_year:

  year to use for out-of-sample validation, overwrites held_out (to
  confirm it does)

- validation_run:

  boolean indicator of whether it's a validation model run or not

- generate_quantities:

  binary vector indicating whether to simulate data from the fitted
  model

  Setting for where to save things

- get_posteriors:

  boolean indicator of whether to return posterior samples

- output_dir:

  output directory, defaults to NULL

- runnumber:

  number to add to runname

- rungroup:

  group to add to runname

  Settings for sampling

- runname:

  name to use for run, if output_dir is NULL

- chains:

  number of chains to run

- iter_sampling:

  number of posterior samples to draw

- iter_warmup:

  number of warmup iterations

- add_sample:

  boolean indicator of whether to return samples

- compile_model:

  boolean indicator of whether to compile the Stan model

- force_recompile:

  boolean indicator of whether to force recompilation of the Stan model

- seed:

  random seed

- refresh:

  number of iterations between progress updates

- adapt_delta:

  target acceptance rate for the No-U-Turn Sampler

- max_treedepth:

  maximum tree depth for the No-U-Turn Sampler

- variational:

  boolean indicator of whether to use variational inference (not yet
  tested)

- nthreads_variational:

  number of threads to use for variational inference

- add_inits:

  boolean indicator of whether to add initial values to the Stan model

- start_year:

  start year of estimates.

- t_star:

  reference year used in model.

- hierarchical_asymptote:

  vector specifying hierarchical structure for asymptote (see Details).

- hierarchical_level:

  vector specifying hierarchical structure for the level in reference
  year (see Details).

- hierarchical_splines:

  vector specifying hierarchical structure for spline coefficients (see
  Details).

- Betas_upper_bound:

  upper bound for the splines parameters

- Betas_lower_bound:

  lower bound for the splines parameters

- Ptilde_low:

  lower bound for the asymptote Ptilde

## Value

fpemplus object.
