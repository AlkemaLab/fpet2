# Plot FP estimates for one country/region, marital status, and indicator

Plot FP estimates for one country/region, marital status, and indicator

## Usage

``` r
plot_estimates_local(
  estimates,
  filtered_data,
  indicator_select,
  estimates2 = NULL,
  estimates3 = NULL,
  estimates4 = NULL,
  modelnames = c("model1", "model2"),
  cols_sourcetypes = c(DHS = "red", DHS0 = "red", MICS = "blue", PMA = "darkgreen", Other
    = "orange", `National survey` = "purple")
)
```

## Arguments

- estimates:

  estimates

- filtered_data:

  observations

- indicator_select:

  indicator

- estimates2:

  optional second model estimates to compare

- estimates3:

  optional third model estimates to compare

- estimates4:

  optional fourth model estimates to compare

- modelnames:

  names of models to use in legend if comparing more than one model

- cols_sourcetypes:

  colors to use for data source types

## Value

ggplot object
