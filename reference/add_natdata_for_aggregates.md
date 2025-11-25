# Add national data for aggregates

This function processes national data and adds it to the aggregate
dataset.

## Usage

``` r
add_natdata_for_aggregates(national_dat_df, regions_dat, is_married)
```

## Arguments

- national_dat_df:

  national survey data prior to processing

- regions_dat:

  region metadata

- is_married:

  boolean

## Value

Processed national data with region_code set to "National"
