# Get Stan data for service statistics data

This function is called from fit_fpem, to prepare service statistics
data for use in model fitting

## Usage

``` r
service_statistics_getstandata(
  service_statistic_df,
  time_index,
  hyper_param,
  geo_unit_index
)
```

## Arguments

- service_statistic_df:

  emus filtered to pop (eg country)

- time_index:

  time_index used in model fit

- hyper_param:

  hyper parameters for EMUs to be used in model fit

- geo_unit_index:

  geo_unit_index used in model fit

## Value

list with dat_emu (for plotting) and emu_list (to pass to stan)
