# Get samples for indicators

This function extracts samples for specified indicators from a posterior
samples object

## Usage

``` r
get_samples_for_indicators(
  indicators,
  has_geo = TRUE,
  samples,
  geo_unit_index = NULL,
  time_index,
  stringtoremove = "",
  rename_modern = FALSE
)
```

## Arguments

- indicators:

  vector of indicator names (e.g., c("contraceptive_use_modern",
  "unmet_need_modern"))

- has_geo:

  logical, if TRUE indicators are \[geo_unit, time\] indicators; if
  FALSE, just \[time\] (for aggregates)

- samples:

  posterior samples object

- geo_unit_index:

  data frame with geo unit information (required if has_geo is TRUE)

- time_index:

  data frame with time information

- stringtoremove:

  string to remove from indicator names

- rename_modern:

  logical, if TRUE renames indicators "modern", "unmet", "trad" to full
  names

## Value

tibble with samples for the specified indicators
