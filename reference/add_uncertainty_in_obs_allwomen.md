# Calculate uncertainty in observations for all-women fit

For all-women FPET fit, adds uncertainty in observations based on
sampling errors and, if specified, also adds NSE uncertainty.

## Usage

``` r
add_uncertainty_in_obs_allwomen(
  samples,
  data_allwomen,
  is_married = TRUE,
  perc_low = 0.025,
  perc_up = 0.975,
  add_nse_uncertainty = TRUE
)
```

## Arguments

- samples:

  posterior::draws object from cmdstanr fit

- data_allwomen:

  list with data for all women used in model fitting

- is_married:

  boolean

- perc_low:

  percentile for lower bound of uncertainty interval, defaults to 0.025

- perc_up:

  percentile for upper bound of uncertainty interval, defaults to 0.975

- add_nse_uncertainty:

  boolean, if TRUE (default), adds NSE uncertainty in addition to
  sampling error uncertainty

## Value

survey data for both marital groups with columns for point estimates and
uncertainty intervals added
