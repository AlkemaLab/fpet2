# Plot FP estimates

Plot FP estimates for multiple countries/regions, marital status, and
indicators

## Usage

``` r
plot_estimates_local_all(
  iso_select = NULL,
  subnational = FALSE,
  results,
  marital_status_all = c("married", "unmarried", "all"),
  indicator_select_all = c("contraceptive_use_modern", "unmet_need_modern", "demand",
    "demand_satisfied_modern", "contraceptive_use_traditional", "contraceptive_use_any",
    "unmet_need_any", "demand_satisfied"),
  dat_emu = NULL,
  results2 = NULL,
  results3 = NULL,
  results4 = NULL,
  modelnames = NULL,
  add_title = TRUE,
  nrow_plot = 2,
  save_plots = TRUE,
  cols_emus = c(clients = "red", visits = "blue", facilities = "darkgreen", users =
    "orange"),
  output_folder = NULL,
  plot_name = "estimates"
)
```

## Arguments

- iso_select:

  name or names of country iso or region code

- subnational:

  logical, whether subnational or national

- results:

  output from fit_fpem()

- marital_status_all:

  choice out of c("married", "unmarried", "all")

- indicator_select_all:

  choice out of c("contraceptive_use_modern", "unmet_need_modern",
  "demand", "demand_satisfied_modern")

- dat_emu:

  data frame of emu data, if NULL, will not plot emu data

- results2:

  optional second model results to compare

- results3:

  optional third model results to compare

- results4:

  optional fourth model results to compare

- modelnames:

  names of models to use in legend if comparing more than one model

- add_title:

  logical, whether to add a title to the plot

- nrow_plot:

  when plotting 4 indicators, use 1 or 2 rows

- save_plots:

  logical, whether to save plots as pdf

- cols_emus:

  colors to use for emu data types

- output_folder:

  folder to save plots in

- plot_name:

  name of pdf file to save plots in

## Value

list of plots for selected countries/regions, marital status, and
indicators
