# Get estimates from global runs

Let's keep this function as we use it for the global runs it worked for
all women but not yet aggregates

## Usage

``` r
get_estimates_globalruns(
  samples,
  geo_unit_index,
  time_index,
  marital_status = "married",
  return_samples = FALSE,
  add_trad = TRUE
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

- marital_status:

  string, "married" or "unmarried"

- return_samples:

  boolean, if TRUE, returns posterior samples along with summaries

- add_trad:

  boolean, if TRUE, includes traditional use in the output

## Value

tibble with estimates (mean and percentiles) for each indicator and
marital status group
