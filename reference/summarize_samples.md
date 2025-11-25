# Summarize samples

Summarize samples

## Usage

``` r
summarize_samples(
  samples,
  percentiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
)
```

## Arguments

- samples:

  tibble with \`value\` column and other grouping columns

- percentiles:

  vector of percentiles to compute

## Value

tibble with summary columns (mean and percentiles) for each group
