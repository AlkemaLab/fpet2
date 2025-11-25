# Use global fit 1a residuals to update assignment of outlier errors

In this function, we calculate residuals from global fit 1a to update
the \`nooutlier\` column for use in further model fitting. In this
function, \`nooutlier\` is set to 0 for observations that are outlying
in 1a, here defined as observations with absolute standardized residuals
above the 90th percentile. By setting \`nooutlier\` to zero, outlier
errors are assigned. (Naming of the column can use an update!)

## Usage

``` r
update_nooutlier_in_global_fit(fit, is_married)
```

## Arguments

- fit:

  global fit 1a

- is_married:

  logical

## Value

data frame with updated nooutlier column.
