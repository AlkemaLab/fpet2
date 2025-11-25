  // data properties, different between marital groups
  int N; // Number of observations
  int S; // Number of sources
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=0, upper=n_geounit> geo_unit; // Geographic unit of each observation, 0 if observation is for an aggregated unit
  // updated to lower of 0 for subnational runs
  array[N] int<lower=0, upper=S> source;   // Source of each observation
  array[N] int<lower=0, upper=1> held_out; // Whether to hold out each observation
  array[N] int<lower=0> isdhs;            // DHS?
  // pma data
  array[N] real<lower=0> ispma;            // PMA?
  array[N] real<lower=0> exact_time;
  int<lower=0> N_pma_country_indices;
  array[N_pma_country_indices] int<lower=1, upper = n_geounit> pma_country_indices;
  array[n_geounit] int npma_c;
  array[n_geounit] int nstart_pma_c;
  int<lower=0> N_pma_obs_indices;
  array[N_pma_obs_indices] int indices_pma;

// for unmarried, add unmarried_
// note that geounit_index need to match between married and unmarried
  int unmarried_N; // Number of observations
  int unmarried_S; // Number of sources
  array[unmarried_N] int<lower=1, upper=T> unmarried_time;     // Time of each observation
  array[unmarried_N] int<lower=0, upper=n_geounit> unmarried_geo_unit; // Geographic unit of each observation, 0 if observation is for an aggregated unit
  // updated to lower of 0 for subnational runs
  array[unmarried_N] int<lower=0, upper=S> unmarried_source;   // Source of each observation
  array[unmarried_N] int<lower=0, upper=1> unmarried_held_out; // Whether to hold out each observation
  array[unmarried_N] int<lower=0> unmarried_isdhs;            // DHS?
  // pma data
  array[unmarried_N] real<lower=0> unmarried_ispma;            // PMA?
  array[unmarried_N] real<lower=0> unmarried_exact_time;
  int<lower=0> unmarried_N_pma_country_indices;
  array[unmarried_N_pma_country_indices] int<lower=1, upper = n_geounit> unmarried_pma_country_indices;
  array[n_geounit] int unmarried_npma_c;
  array[n_geounit] int unmarried_nstart_pma_c;
  int<lower=0> unmarried_N_pma_obs_indices;
  array[unmarried_N_pma_obs_indices] int unmarried_indices_pma;

  matrix[n_geounit, T] prop_married_rt; // proportion married for each region-year
  vector[T] prop_married_t; // at national level, same as prop_married_rt for national run

