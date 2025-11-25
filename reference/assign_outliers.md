# Assign outliers

Function to assign outliers based on a set of rules. If columns
possible_outlier or possible_outlier_userinput are present, these are
used to define outliers.

## Usage

``` r
assign_outliers(fp_dat, outlier_record_ids = NULL, married)
```

## Arguments

- fp_dat:

  survey data

- outlier_record_ids:

  record_ids of outliers (if any)

- married:

  boolean, needed when bias needs to be assessed

## Value

data frame with new column nooutlier (1 = not an outlier, 0 = outlier)
