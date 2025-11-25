FPET2
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 

R package for fitting Bayesian transition models to produce estimates
and projections of family planning indicators. This work was supported,
in whole or in part, by the Bill & Melinda Gates Foundation (INV-00844).

# Installation

Dependencies

- `cmdstanr`: Instructions for installing `cmdstanr` are available in
  their [Getting
  started](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) guide.
- R package `localhierarchy`, available at
  [github.com/AlkemaLab/localhierarchy](github.com/AlkemaLab/localhierarchy).
  You can install it
  using`remotes::install_github("AlkemaLab/localhierarchy")`.

Install `fpet2` from Github:

`remotes::install_github("AlkemaLab/fpet2")`

# Information

The modeling approach is described in Alkema et al. (2024). This
approach builds off collaborative research with the United Nations
Population Division (UNPD). References to earlier and related work are
included below.

The national-level data and population estimates by marital status are
compiled and published by the United Nations Population Division (UNPD).
Track20 compiled a database of subnational data.

References that describe the current approach:

- FPET overview: L. Alkema, H. Susmann, E. Ray, S. Mooney, N. Cahill, K.
  Bietsch, A. Jayachandran, et al. (2024). Statistical Demography Meets
  Ministry of Health: The Case of the Family Planning Estimation Tool.
  See <https://doi.org/10.48550/arXiv.2501.00007>.

- Use of service statistics data (Estimated modern use): S. Mooney, L.
  Alkema, E. Sonneveldt, K. Bietsch, J. Williamson, N. Cahill (2024).
  Enhancing the Use of Family Planning Service Statistics Using a
  Bayesian Modelling Approach to Inform Estimates of Modern
  Contraceptive Use in Low- and Middle-Income Countries. See
  <https://doi.org/10.48550/arXiv.2412.08606>.

References that describe earlier approaches:

- Alkema, L. et al. (2013). National, regional and global rates and
  trends in contraceptive prevalence and unmet need for family planning
  between 1990 and 2015: a systematic and comprehensive analysis. The
  Lancet. Available at
  <http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(12)62204-1/abstract>.

- Cahill et al., (2017). Modern contraceptive use, unmet need, and
  demand satisfied among women of reproductive age who are married or in
  a union in the focus countries of the Family Planning 2020 initiative:
  a systematic analysis using the Family Planning Estimation Tool. The
  Lancet. Available at
  <http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(17)33104-5/fulltext>.

- Kantorova et al. (2020). Estimating progress towardsmeeting women’s
  contraceptive needs in 185 countries: A Bayesian hierarchical
  modellingstudy. PLOS Medicine, 17(2):e1003026. Available at
  <https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1003026>.

- Wheldon et al. (2018). Methods for estimating and projecting key
  family planning indicators among all women of reproductive age. United
  Nations, Department of Economic and Social Affairs, Population
  Division, Technical Paper No. 2. New York: United Nations. Available
  at
  <https://www.un.org/en/development/desa/population/publications/pdf/technical/TP2018-2.pdf>.

The model is written in R and Stan. For any questions, feedback or
suggestions, please contact Leontine at lalkema(at)umass.edu.

# Example

See articles on package website, at
<a href="https://alkemalab.github.io/fpet2/index.html"
class="uri">https://alkemalab.github.io/fpet2</a>.

# Citation

Please cite as follows:

``` r
citation("fpet2")
#> To cite package 'fpet2' in publications use:
#> 
#>   Alkema L, Susmann H, Mooney S, Ray E (2025). _fpet2: An R Package for
#>   Estimation and Forecasting of Family Planning Indicators_. R package
#>   version 1.2, <https://github.com/AlkemaLab/fpet2>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {fpet2: An R Package for Estimation and Forecasting of Family Planning Indicators},
#>     author = {Leontine Alkema and Herbert Susmann and Shauna Mooney and Evan Ray},
#>     year = {2025},
#>     note = {R package version 1.2},
#>     url = {https://github.com/AlkemaLab/fpet2},
#>   }
```
