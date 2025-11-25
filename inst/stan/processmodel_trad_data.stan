
  int<lower=0, upper=1> z_fix_smoothing;

  int<lower=0> z_Omega_raw_n_terms;
  int<lower=0> z_Omega_n_sigma;
  int<lower=0> z_Omega_raw_n_terms_fixed;
  int<lower=0> z_Omega_raw_n_terms_estimate;
  int<lower=0> z_Omega_n_sigma_fixed;
  int<lower=0> z_Omega_n_sigma_estimate;
  array[z_Omega_n_sigma + 1] int<lower=1, upper=z_Omega_raw_n_terms> z_Omega_re_start;
  array[z_Omega_n_sigma + 1] int<lower=1, upper=z_Omega_raw_n_terms> z_Omega_re_end;
  matrix[n_geounit, z_Omega_raw_n_terms] z_Omega_model_matrix;
  real<lower = verysmallnumber> z_Omega_scalarprior_sd;
  real z_Omega_scalarprior_mean;
  real<lower = verysmallnumber> z_Omega_prior_sd_sigma_estimate;
  vector[z_Omega_raw_n_terms_fixed] z_Omega_raw_fixed;
  vector<lower=0>[z_Omega_n_sigma_fixed] z_Omega_sigma_fixed;


  // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[z_fix_smoothing] real<lower=0, upper=1> z_Rho_fixed;
  array[z_fix_smoothing] real<lower=0> z_Tau_fixed;

