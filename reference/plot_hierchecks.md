# Plot hierarchical checks

Plot hierarchical checks

## Usage

``` r
plot_hierchecks(
  fit,
  is_married,
  add_priors_transitions = TRUE,
  add_priors_trad = FALSE,
  add_trad = FALSE
)
```

## Arguments

- fit:

  fit object from a global fit

- is_married:

  logical, is this a married model?

- add_priors_transitions:

  Plot prior-posts for sigmas?

- add_priors_trad:

  Plot prior-posts for trad sigmas?

- add_trad:

  Plot trad parameters?

## Value

NULL, but saves a pdf with plots to the output directory
