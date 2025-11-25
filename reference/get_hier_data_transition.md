# Get hierarchical data for transition model

Get hierarchical data for transition model

## Usage

``` r
get_hier_data_transition(
  geo_unit_index,
  prefix = "",
  stan_data_settings,
  stan_spline_data,
  hierarchical_terms_and_fixed,
  global_fit = NULL
)
```

## Arguments

- geo_unit_index:

  geo_unit_index for model fit

- prefix:

  prefix to add to parameter names, e.g., "d\_" or "z\_"

- stan_data_settings:

  list of stan data settings, e.g., bounds and priors

- stan_spline_data:

  stan data for splines, needed for dimension of splines coeff vector

- hierarchical_terms_and_fixed:

  list with info on levels and what's fixed from other function

- global_fit:

  fit object from a global fit, or NULL if not doing a local fit

## Value

list with hier_stan_data_settings, hier_stan_data, and hier_data
