# Fit FPET national or subnational local models for ALL women, with or without EMUs.

Fit FPET national or subnational local models for ALL women, with or
without EMUs.

## Usage

``` r
fit_fpem(
  survey_df,
  population_df,
  service_statistic_df = NULL,
  subnational = FALSE,
  national_dat_df = NULL,
  chains = 4,
  iter_sampling = 200,
  iter_warmup = 150,
  seed = 1234,
  regions_dat,
  end_year = 2030
)
```

## Arguments

- survey_df:

  Survey data for model fitting, filtered to pop of interest

- population_df:

  Population data, filtered to pop of interest

- service_statistic_df:

  EMU data, filtered to pop of interest

- subnational:

  logical, TRUE if fitting to subnational data. For subnational fits,
  survey_df may contain all or just 1 region.

- national_dat_df:

  National data to be used/displayed for subnational runs, filtered to
  pop of interest. Model is fitted to surveys in national_dat_df that
  are not in survey_df.

- chains:

  Number of chains for Stan

- iter_sampling:

  Number of iterations for sampling

- iter_warmup:

  Number of iterations for warmup

- seed:

  Random seed for Stan

- regions_dat:

  Data frame with region information

- end_year:

  End year for the model, default is 2030

## Value

results list with estimates, samples, data_allwomen, observations, and
dat_emu = dat_emu
