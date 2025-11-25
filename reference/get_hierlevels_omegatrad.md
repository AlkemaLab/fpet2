# Get list with information on hierarchical levels and which ones to fix for omegatrad

Get list with information on hierarchical levels and which ones to fix
for omegatrad

## Usage

``` r
get_hierlevels_omegatrad(
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

  character, prefix for process model (should be "z\_" for traditional)

- add_subnational_hierarchy:

  character, name of subnational hierarchy level to add, e.g., "subnat"

## Value

list with hierarchical terms and what's fixed
