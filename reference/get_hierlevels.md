# Get list with information on hierarchical levels for transition model parameters, and which ones to fix

Get list with information on hierarchical levels for transition model
parameters, and which ones to fix

## Usage

``` r
get_hierlevels(
  runstep,
  hierarchical_terms = NULL,
  global_fit = NULL,
  prefix = "",
  add_subnational_hierarchy = NULL
)
```

## Arguments

- runstep:

  character, from model fit

- hierarchical_terms:

  list of hierarchical terms from fit_model inputs or NULL

- global_fit:

  fit object from a global fit

- prefix:

  character, prefix for process model, e.g., "", d\_", "z\_"

- add_subnational_hierarchy:

  character, name of subnational hierarchy level to add, e.g., "subnat"

## Value

list with hierarchical terms and what's fixed
