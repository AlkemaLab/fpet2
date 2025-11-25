# Get stan data for non-sampling errors (nonse)

Get stan data for non-sampling errors (nonse)

## Usage

``` r
get_nonse_standata(fix_nonse, global_fit, source_index, prefix = "DM1_")
```

## Arguments

- fix_nonse:

  logical, whether to fix nonse parameters to values from a global fit

- global_fit:

  fit object from a global fit, if fix_nonse is TRUE

- source_index:

  source_index data frame for the local fit

- prefix:

  prefix for nonse parameters, defaults to "DM1\_"

## Value

list of stan data elements for nonse
