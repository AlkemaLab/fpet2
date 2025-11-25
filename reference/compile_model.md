# Compile a Stan model

Compile a Stan model

## Usage

``` r
compile_model(
  variational,
  nthreads_variational,
  force_recompile,
  stan_file_path
)
```

## Arguments

- variational:

  boolean

- nthreads_variational:

  if variational, number of threads to use

- force_recompile:

  boolean

- stan_file_path:

  path to the stan file

## Value

cmdstanr::cmdstan_model object
