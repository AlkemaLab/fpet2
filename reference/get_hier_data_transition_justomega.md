# Get hierarchical data for Omegas in transition model

Perhaps to combine with \`get_hier_data_transition\`

## Usage

``` r
get_hier_data_transition_justomega(
  geo_unit_index,
  prefix = "",
  stan_data_settings,
  hierarchical_terms_and_fixed,
  global_fit = NULL,
  runstep
)
```

## Arguments

- geo_unit_index:

  geo_unit_index for model fit

- prefix:

  prefix to add to parameter names, e.g., "d\_" or "z\_"

- stan_data_settings:

  list of stan data settings, e.g., bounds and priors

- hierarchical_terms_and_fixed:

  list with info on levels and what's fixed from other function

- global_fit:

  fit object from a global fit, or NULL if not doing a local fit

- runstep:

  what model is fitted

## Value

list with hier_stan_data_settings, hier_stan_data, and hier_data
