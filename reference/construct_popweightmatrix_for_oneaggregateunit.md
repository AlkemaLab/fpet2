# Construct population weight matrix for one aggregate unit

Constructs a row vector of population weights for regions in a single
aggregate unit.

## Usage

``` r
construct_popweightmatrix_for_oneaggregateunit(
  aggregation_meta_one_unit,
  population_selectedyear,
  name_region_index = "region_index"
)
```

## Arguments

- aggregation_meta_one_unit:

  tibble with columns agg_unit_name and region_index

- population_selectedyear:

  tibble with region_index and population, needs to include all regions

## Value

row vector with weights for each region in the aggregate unit (0 if
region is not included)
