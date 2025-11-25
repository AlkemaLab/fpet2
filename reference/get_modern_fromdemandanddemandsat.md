# Get modern from demand and demand satisfied

Model returns Eta (demandsat) and d_Eta (demand) \# real demand =
inv_tr_eta(tr_d_Eta_obs\[geo_unit\[i\], time\[i\]\]); \# modern\[i\] =
inv_tr_eta(tr_Eta_obs\[geo_unit\[i\], time\[i\]\]) \* demand; \#
logit_unmetovernonmodern\[i\] = logit((demand -
modern\[i\])/(1-modern\[i\]));

## Usage

``` r
get_modern_fromdemandanddemandsat(demand, demand_satisfied)
```

## Arguments

- demand:

  demand

- demand_satisfied:

  demand satisfied (modern/demand)

## Value

modern
