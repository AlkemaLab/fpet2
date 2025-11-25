# Get standard errors for FP indicators in sequential DM approach

Get standard errors for FP indicators in sequential DM approach

## Usage

``` r
get_se_sequentialdm(nat_data, use_small_SE_setting = TRUE)
```

## Arguments

- nat_data:

  data frame with survey data

- use_small_SE_setting:

  logical, if TRUE (default), use smaller effective sample sizes when
  imputing SEs

## Value

tibble with SEs added
