# Get estimates for all-women model fits

This function processes posterior samples from FPET to get estimates for
different marital status groups and indicators. It can also handle
subnational estimation if specified.

## Usage

``` r
get_estimates_fpem(
  samples,
  geo_unit_index,
  time_index,
  subnational = FALSE,
  return_samples = FALSE
)
```

## Arguments

- samples:

  tibble of posterior samples with columns for indicators,
  geo_unit_index, time_index, and value

- geo_unit_index:

  geo_unit_index from model fit

- time_index:

  time_index from model fit

- subnational:

  boolean, if TRUE, also computes national aggregates for subnational
  models

- return_samples:

  boolean, if TRUE, returns posterior samples along with summaries

## Value

tibble with estimates (mean and percentiles) for each indicator and
marital status group
