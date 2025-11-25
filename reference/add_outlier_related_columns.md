# Add possible-outlier related columns

To survey data, add columns possible_outliers and
possible_outlier_userinput. Possible outliers are assessed using the
\`assign_outliers\` function.

## Usage

``` r
add_outlier_related_columns(dat)
```

## Arguments

- dat:

  survey data

## Value

data frame with new columns possible_outlier (1 = possible outlier, 0 =
not a possible outlier) and possible_outlier_userinput (NA)
