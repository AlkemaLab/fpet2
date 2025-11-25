# Get smoothing standata

Get smoothing standata

## Usage

``` r
get_smoothing_standata(fix_smoothing, global_fit, prefix = "")
```

## Arguments

- fix_smoothing:

  logical, whether to fix smoothing parameters

- global_fit:

  fit object from a global fit, used to extract smoothing parameters if
  fix_smoothing is TRUE

- prefix:

  prefix to add to the parameter names, defaults to ""

## Value

list with smoothing standata
