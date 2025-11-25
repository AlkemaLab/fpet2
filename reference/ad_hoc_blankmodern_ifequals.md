# Blank out modern use if modern cp and any CP are equal

Blank out modern use if modern cp and any CP are equal

## Usage

``` r
ad_hoc_blankmodern_ifequals(obs)
```

## Arguments

- obs:

  data frame with contraceptive_use_any, contraceptive_use_modern,
  contraceptive_use_traditional

## Value

data frame with contraceptive_use_modern set to NA if it equals
contraceptive_use_any and contraceptive_use_traditional is NA
