# probit \<- function(x) pnorm(x) \# stan \#phiapprox \<- function(x) inv_logit(0.07056 \* x ^ 3 + 1.5976 \* x) inv_probit \<- function(x) qnorm(x)

\#' Get standard error on inverse probit scale from proportion and its
standard error \#' \#' @param prop proportion \#' @param se_prop
standard error of the proportion (original scale) \#' \#' @returns
standard error on inverse probit scale \#' @keywords internal
get_se_invprobitprop \<- function(prop, se_prop) \#
se_prop\*abs(derivate evaluated at prop)
se_prop\*1/dnorm(inv_probit(prop)) Logit function

## Usage

``` r
logit(x)
```

## Arguments

- x:

  input value

## Value

logit of x
