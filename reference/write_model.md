# Write stan models

Write stan models

## Usage

``` r
write_model(
  add_aggregates = FALSE,
  add_trad = TRUE,
  all_women = FALSE,
  add_emu = FALSE
)
```

## Arguments

- add_aggregates:

  Add national aggregates? makes sense only for subnational all-region
  run for 1 country

- add_trad:

  Include trad use estimation and data?

- all_women:

  Estimate for all_women or by marital group?

- add_emu:

  Include EMU estimation and data?

## Value

writes a stan file to the stan directory. The model files are named
based on the options chosen.
