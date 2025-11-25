# Workflow for FPET global fitting for national-level estimation

## Introduction

In this vignette we illustrate how to do global fitting for national
estimation of family planning indicators using the FPET2 R package.

``` r
library(fpet2)
options(cmdstanr_warn_inits = FALSE)
```

## Global survey data base

For global fitting, we use a global data base of survey data. We also
use meta information on countries and world regions.

Read data

``` r
data_folder <- here::here("data-raw")
survey_data_file <- file.path(data_folder, "Track20 2025 Database for FPET2 051325.csv")
survey_df_all <- readr::read_csv(survey_data_file, show_col_types = FALSE)
```

Global national model fitting is done separately per marital group. We
illustrate it here for married women.

We process the data for that marital group:

``` r
is_married <- TRUE
dat <- process_data(survey_df_all,
                    regions_dat = fpet2::regions_all,
                    is_married = is_married)
```

    ## [1] "When imputing SEs for DHS, we use the effective sample size of its preceding survey"

    ## Joining with `by = join_by(record_id_fixed)`
    ## Joining with `by = join_by(record_id_fixed)`

    ## [1] "We define possible outliers based on column possible_outlier"

## Model settings: MCMC settings and output directory

For model fitting, we define MCMC settings. We can choose a test setting
or a longer run. We use a test setting here for illustration only. NOTE
THAT THIS WILL PRODUCE WARNINGS AND RESULTS THAT ARE NOT TO BE USED!!

``` r
# test settings
chains = 4
iter_sampling = 5
iter_warmup = 5

# longer run
# chains = 16
# iter_sampling = 50
# iter_warmup = 150
```

When fitting the model, outputs are saved in an output directory. There
are different options to define this directory. By default, the
directory is automatically created one directory up from this package’s
directory, in a directory called `bayestransition_output`. Alternatively
(and easier for this article), the user can specify an output directory
to use:

``` r
# Create an output directory in a temporary location
output_dir <- file.path(tempdir(), "vignette_output")
dir.create(output_dir, showWarnings = FALSE)

# Print the path so the user knows where to find it
cat("Outputs are saved in:", output_dir, "\n")
```

    ## Outputs are saved in: /tmp/Rtmp8bOERs/vignette_output

## Model fitting

The `fit_model` function is used for global model fitting. Its help
function explains its use. Its argument `run_step` is used to define the
type of model fitting. For global national fitting, we do the following
fits:

- `runstep = "step1a"`: Get longterm trends by fitting a model witout
  AR(1) deviation terms;

- `runstep = "step1b"`: Estimate parameters governing country-specific
  variations in transition functions, short-term fluctuations and data
  quality;

- `runstep = "step2"`: Estimate parameters associated with traditional
  use.

Let’s start with 1a:

``` r
fit1a <- fit_model(runstep = "step1a",
                   survey_df = dat,
                   is_married = is_married,
                   get_posteriors = FALSE,
                   output_dir = output_dir, 
                   chains = chains,
                   iter_sampling =  iter_sampling,
                   iter_warmup = iter_warmup)
```

    ## [1] "For 1a, we take all levels from the fit_model inputs"
    ## [1] "We do not fix any terms or sigmas of hierarchical models."
    ## [1] "For 1a, we take all levels from the fit_model inputs"
    ## [1] "We do not fix any terms or sigmas of hierarchical models."
    ## [1] "We do give nonSE to DHS (by temporarily renaming DHS into DHS0)"
    ## [1] "output directory is /tmp/Rtmp8bOERs/vignette_output"
    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 WARNING: No variance estimation is 
    ## Chain 1          performed for num_warmup < 20 
    ## Chain 2 WARNING: No variance estimation is 
    ## Chain 2          performed for num_warmup < 20 
    ## Chain 3 WARNING: No variance estimation is 
    ## Chain 3          performed for num_warmup < 20 
    ## Chain 4 WARNING: No variance estimation is 
    ## Chain 4          performed for num_warmup < 20 
    ## Chain 3 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 4 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 1 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 2 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 2 finished in 0.5 seconds.
    ## Chain 1 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 1 finished in 0.8 seconds.
    ## Chain 3 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 3 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 3 finished in 3.2 seconds.
    ## Chain 4 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 4 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 4 finished in 13.4 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 4.5 seconds.
    ## Total execution time: 14.0 seconds.

    ## Warning: 20 of 20 (100.0%) transitions ended with a divergence.
    ## See https://mc-stan.org/misc/warnings for details.

Outputs are stored in

``` r
fit1a$output_dir
```

    ## [1] "/tmp/Rtmp8bOERs/vignette_output"

The next fit 1b needs summary information from 1a passed on through its
`global_fit` object. Here we pass the one we just fitted. For fit1b, we
also update the outlier classifications in the survey data used, to
allow for data points that are outlying in fit 1a to be possibly
outlying, as follows:

``` r
dat_updated <- update_nooutlier_in_global_fit(fit1a, is_married)
```

    ## Warning: Dropping 'draws_df' class as required metadata was removed.
    ## Warning: Dropping 'draws_df' class as required metadata was removed.

    ## Joining with `by = join_by(i)`
    ## Joining with `by = join_by(c)`

Now let’s fit step 1b:

``` r
fit1b <- fit_model(runstep = "step1b",
                   global_fit = fit1a,
                   survey_df = dat_updated,
                   is_married = is_married,
                    output_dir = output_dir, 
                   get_posteriors = FALSE,
                   chains = chains,
                   iter_sampling =  iter_sampling,
                   iter_warmup = iter_warmup
)
```

    ## [1] "We use a global fit, and take selected settings from there."
    ## [1] "settings for the spline_degree and num_knots taken from global run"
    ## [1] "Setting for tstar taken from global run"
    ## [1] "We take all hier terms from the global fit, using prefix"
    ## [1] "For hierarchical terms, we fix things up to the 2nd-lowest level."
    ## [1] "For sigma terms in hierarchical models for demand and ds, we fix things up to the 2nd-lowest level."
    ## [1] "We take all hier terms from the global fit, using prefix"
    ## [1] "For hierarchical terms, we fix things up to the 2nd-lowest level."
    ## [1] "For sigma terms in hierarchical models for demand and ds, we fix things up to the 2nd-lowest level."
    ## [1] "output directory is /tmp/Rtmp8bOERs/vignette_output"
    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 WARNING: No variance estimation is 
    ## Chain 1          performed for num_warmup < 20 
    ## Chain 2 WARNING: No variance estimation is 
    ## Chain 2          performed for num_warmup < 20 
    ## Chain 3 WARNING: No variance estimation is 
    ## Chain 3          performed for num_warmup < 20 
    ## Chain 4 WARNING: No variance estimation is 
    ## Chain 4          performed for num_warmup < 20 
    ## Chain 1 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 3 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 4 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 2 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 2 finished in 1.9 seconds.
    ## Chain 3 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 3 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 3 finished in 3.2 seconds.
    ## Chain 4 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 1 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 4 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 4 finished in 3.7 seconds.
    ## Chain 1 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 1 finished in 4.2 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 3.3 seconds.
    ## Total execution time: 4.7 seconds.

    ## Warning: 20 of 20 (100.0%) transitions ended with a divergence.
    ## See https://mc-stan.org/misc/warnings for details.

In step 2, the parameters for traditional use are estimated:

``` r
fit2 <- fit_model(runstep = "step2",
                   global_fit = fit1b,
                   survey_df = dat_updated,
                   is_married = is_married,
                   get_posteriors = FALSE,
                  output_dir = output_dir, 
                  chains = chains,
                  iter_sampling =  iter_sampling,
                  iter_warmup = iter_warmup)
```

    ## [1] "We use a global fit, and take selected settings from there."
    ## [1] "settings for the spline_degree and num_knots taken from global run"
    ## [1] "Setting for tstar taken from global run"
    ## [1] "We take all hier terms from the global fit, using prefix"
    ## [1] "For hierarchical terms, we fix things up to the 2nd-lowest level."
    ## [1] "We fix all sigmas of hierarchical models for demand and ds."
    ## [1] "We take all hier terms from the global fit, using prefix"
    ## [1] "For hierarchical terms, we fix things up to the 2nd-lowest level."
    ## [1] "We fix all sigmas of hierarchical models for demand and ds."
    ## [1] "For 2, we take all levels for trad from the fit_model inputs"
    ## [1] "We do not fix any terms or sigmas of hierarchical models."
    ## [1] "output directory is /tmp/Rtmp8bOERs/vignette_output"
    ## Running MCMC with 4 parallel chains...
    ## 
    ## Chain 1 WARNING: No variance estimation is 
    ## Chain 1          performed for num_warmup < 20 
    ## Chain 2 WARNING: No variance estimation is 
    ## Chain 2          performed for num_warmup < 20 
    ## Chain 3 WARNING: No variance estimation is 
    ## Chain 3          performed for num_warmup < 20 
    ## Chain 4 WARNING: No variance estimation is 
    ## Chain 4          performed for num_warmup < 20 
    ## Chain 1 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 3 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 4 Iteration: 1 / 10 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 1 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 4 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 1 finished in 1.5 seconds.
    ## Chain 4 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 4 finished in 1.9 seconds.
    ## Chain 2 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 2 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 2 finished in 4.8 seconds.
    ## Chain 3 Iteration: 6 / 10 [ 60%]  (Sampling) 
    ## Chain 3 Iteration: 10 / 10 [100%]  (Sampling) 
    ## Chain 3 finished in 6.3 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 3.6 seconds.
    ## Total execution time: 7.1 seconds.

    ## Warning: 20 of 20 (100.0%) transitions ended with a divergence.
    ## See https://mc-stan.org/misc/warnings for details.

We can plot all estimates using the `plot_estimates_local` function.
Relevant code is added below. We are not showing results hereto avoid
confusion, given that test settings were used, so the results are
incorrect.

Let’s create the results first:

``` r
results2 <- list()
results2$estimates <- get_estimates_globalruns(
  samples = fit2$samples, 
  geo_unit_index = fit2$geo_unit_index,
  time_index = fit2$time_index,
  marital_status = ifelse(is_married, "married", "unmarried"),
  return_samples = FALSE)
```

Add uncertainty in observations (so that we can display the 95% CIs
associated with each data point):

``` r
results2$observations <- add_uncertainty_in_obs(fit  = fit2)
```

Now let’s plot! (again, not showing the results…)

``` r
plots <-
  plot_estimates_local_all(
      results = results2,
      marital_status_all =  ifelse(is_married, "married", "unmarried"),
      indicator_select_all =  c("contraceptive_use_modern",
                 "unmet_need_modern",
                 "demand", #"demand_satisfied_modern",
                 "contraceptive_use_traditional"),
      save_plots = FALSE)
# plots are in a list, eg check 
plots[[1]][[1]]
```

## Local fitting

For local fitting, the `fit_model` function can be used as well, using
`runstep = "local_national"`. This results in local fitting for one
marital group. For local fitting for all women, the function `fit_fpem`
can be used (see other article).
