# Fit model to produce all women estimates

This function is called from fit_fpem, to do the model fitting for all
women

## Usage

``` r
fit_model_allwomen(
  data_allwomen,
  runstep = "local_national",
  population_df,
  emu_list = NULL,
  end_year,
  chains = 4,
  iter_sampling = 200,
  iter_warmup = 150,
  add_inits = FALSE,
  add_sample = TRUE,
  compile_model = TRUE,
  force_recompile = FALSE,
  seed = 1234,
  refresh = 50,
  adapt_delta = 0.9,
  max_treedepth = 14,
  variational = FALSE,
  nthreads_variational = 8
)
```

## Arguments

- data_allwomen:

  list with input data for married and data for unmarried women

- runstep:

  local_national or local_subnational

- population_df:

  pop counts for each geounit-year combi

- emu_list:

  emu info from \`service_statistics_getstandata\`

- chains:

  Number of chains for Stan

- iter_sampling:

  Number of iterations for sampling

- iter_warmup:

  Number of iterations for warmup

- add_inits:

  Use function to sample initial values? Defaults to FALSE

- add_sample:

  TRUE

- compile_model:

  TRUE

- force_recompile:

  FALSE

- seed:

  Random seed for Stan

- refresh:

  10

- adapt_delta:

  0.9

- max_treedepth:

  14

- variational:

  FALSE

- nthreads_variational:

  8

## Value

results list with stan_model, stan_data, samples, runstep
