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
