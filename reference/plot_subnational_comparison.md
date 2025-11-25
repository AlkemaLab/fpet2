# Plot subnational comparison

Plot subnational comparison

## Usage

``` r
plot_subnational_comparison(
  results,
  indicator_select,
  year_select,
  marital_status_select,
  ymin_select = 0,
  ymax_select = NA,
  arrange_point = TRUE
)
```

## Arguments

- results:

  output from fit_fpem()

- indicator_select:

  choice out of c("contraceptive_use_modern", "unmet_need_modern",
  "demand", "demand_satisfied_modern")

- year_select:

  year to plot

- marital_status_select:

  choice out of c("all", "married", "unmarried")

- ymin_select:

  default 0

- ymax_select:

  default NA

- arrange_point:

  logical, if TRUE arrange points by value of indicator_select

## Value

a ggplot object
