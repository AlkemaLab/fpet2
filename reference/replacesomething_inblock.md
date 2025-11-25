# Replace something in block

Replace something in block

## Usage

``` r
replacesomething_inblock(
  block,
  params = c("Epsilon", "Eta", "Betas", "Ptilde", "Omega", "Rho", "Tau"),
  prefix = "d_"
)
```

## Arguments

- block:

  string with stan code block

- params:

  vector of parameter names to replace

- prefix:

  prefix to add to parameter names

## Value

modified block
