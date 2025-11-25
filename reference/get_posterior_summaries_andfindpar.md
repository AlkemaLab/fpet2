# Get posterior summaries for parameters to fix

This function collects posterior summaries for a set of parameters

## Usage

``` r
get_posterior_summaries_andfindpar(
  fit,
  process_indicator_prefixes = c("", "d_", "z_"),
  dm_indicator_prefixes = c("DM1_", "DM2_", "DM3_"),
  process_params = c("Betas_raw", "Betas_sigma", "Ptilde_raw", "Ptilde_sigma",
    "Omega_raw", "Omega_sigma"),
  process_smoothing_params = c("Rho", "Tau"),
  process_subnat_smoothing_params = c("rho_correlationeps"),
  dm_params = c("nonse", "sdbias", "rho_pma"),
  dm_outlier_params = c("global_shrinkage_dm", "caux_dm")
)
```

## Arguments

- fit:

  fit object from a global fit

- process_indicator_prefixes:

  defaults to c("", "d\_", "z\_")

- dm_indicator_prefixes:

  defaults to c("DM1\_", "DM2\_", "DM3\_")

- process_params:

  defaults to \_raw and \_sigma for Betas, Omega, and Ptilde

- process_smoothing_params:

  defaults to c("Rho", "Tau")

- process_subnat_smoothing_params:

  defaults to c("rho_correlationeps")

- dm_params:

  defaults to c("nonse", "sdbias", "rho_pma")

- dm_outlier_params:

  defaults to c("global_shrinkage_dm", "caux_dm")

## Value

tibble with posterior summaries for parameters to fix
