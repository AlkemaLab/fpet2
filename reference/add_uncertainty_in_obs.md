# Calculate uncertainty in observations for fit to one marital group

For one marital group, adds uncertainty in observations based on
sampling errors and, if specified, also adds NSE uncertainty. To be used
for global fits for one marital group.

## Usage

``` r
add_uncertainty_in_obs(
  fit,
  samples_scale = NULL,
  perc_low = 0.025,
  perc_up = 0.975,
  add_nse_uncertainty = TRUE
)
```

## Arguments

- fit:

  list with fit information, including data and posterior samples

- samples_scale:

  if provided, used instead of fit\$samples\$draws("DMX_scale")

- perc_low:

  percentile for lower bound of uncertainty interval, defaults to 0.025

- perc_up:

  percentile for upper bound of uncertainty interval, defaults to 0.975

- add_nse_uncertainty:

  boolean, if TRUE (default), adds NSE uncertainty in addition to
  sampling error uncertainty

- is_married:

  boolean

## Value

survey data with columns for point estimates and uncertainty intervals
added
