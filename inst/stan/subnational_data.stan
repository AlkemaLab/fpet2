  // subnational
  int n_agg_units;                            // Number of geo units observed at aggregated spatial scale (e.g. number of nations)
  array[N] int<lower=0, upper=1> is_agg_obs; // Whether each observation is for a lowest-level geo unit or aggregated geo units
  array[N] int<lower=0, upper=n_agg_units> agg_unit; // For observations at aggregated spatial scales, index of aggregated unit, 0 if observation is not for an aggregated unit
  array[n_agg_units, n_geounit, T] real geo_unit_pop_wt; // Weights used to aggregate geo units. Each column should sum to 1
  // subnational indicator determines whether aggregates are calculated or not
  int<lower=0, upper=1> subnational;
  // 0: smoothing terms independent across lowest-level geo units
  // 1: smoothing terms correlated among lowest-level geo units in the same group
  int<lower=0, upper=1> correlated_smoothing;
  int n_cor_smoothing_blocks; // basically, number of national level units
  // number of low-level geo units in each block
  array[n_cor_smoothing_blocks] int<lower=1, upper=n_geounit> cor_smoothing_block_sizes;
  int max_cor_smoothing_block_size;

  int<lower=0, upper=1> fix_subnat_corr;
  array[fix_subnat_corr ? 1 : 0]  real<lower=0, upper = 1> rho_correlationeps_fixed;
