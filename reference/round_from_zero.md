# Round from zero

Taken from UNPD/fpemlocal:
https://github.com/AlkemaLab/fpemlocal/blob/master/R/rounding_and_recalculation.R

## Usage

``` r
round_from_zero(x)
```

## Arguments

- x:

  proportion to round

## Value

proportions less than 0.001 are rounded upwards to 0.001
