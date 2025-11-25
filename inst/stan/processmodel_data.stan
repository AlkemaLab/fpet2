  real Betas_lower_bound;
  real Betas_upper_bound;
  real Ptilde_low;

  int<lower=0> Ptilde_raw_n_terms;
  int<lower=0> Ptilde_raw_n_terms_fixed;
  int<lower=0> Ptilde_raw_n_terms_estimate;
  int<lower=0> Ptilde_n_sigma;
  int<lower=0> Ptilde_n_sigma_fixed;
  int<lower=0> Ptilde_n_sigma_estimate;
  array[Ptilde_n_sigma + 1] int<lower=1, upper=Ptilde_raw_n_terms> Ptilde_re_start;
  array[Ptilde_n_sigma + 1] int<lower=1, upper=Ptilde_raw_n_terms> Ptilde_re_end;
  matrix[n_geounit, Ptilde_raw_n_terms] Ptilde_model_matrix;
  real<lower = verysmallnumber> Ptilde_scalarprior_sd;
  real Ptilde_scalarprior_mean;
  real<lower = verysmallnumber> Ptilde_prior_sd_sigma_estimate;
  vector[Ptilde_raw_n_terms_fixed] Ptilde_raw_fixed;
  vector<lower=0>[Ptilde_n_sigma_fixed] Ptilde_sigma_fixed;

  int<lower=0> Omega_raw_n_terms;
  int<lower=0> Omega_n_sigma;
  int<lower=0> Omega_raw_n_terms_fixed;
  int<lower=0> Omega_raw_n_terms_estimate;
  int<lower=0> Omega_n_sigma_fixed;
  int<lower=0> Omega_n_sigma_estimate;
  array[Omega_n_sigma + 1] int<lower=1, upper=Omega_raw_n_terms> Omega_re_start;
  array[Omega_n_sigma + 1] int<lower=1, upper=Omega_raw_n_terms> Omega_re_end;
  matrix[n_geounit, Omega_raw_n_terms] Omega_model_matrix;
  real<lower = verysmallnumber> Omega_scalarprior_sd;
  real Omega_scalarprior_mean;
  real<lower = verysmallnumber> Omega_prior_sd_sigma_estimate;
  vector[Omega_raw_n_terms_fixed] Omega_raw_fixed;
  vector<lower=0>[Omega_n_sigma_fixed] Omega_sigma_fixed;

  int Betas_k_terms; // same as k; to do clean to have just 1
  int<lower=0> Betas_raw_n_terms;
  int<lower=0> Betas_raw_n_terms_fixed;
  int<lower=0> Betas_raw_n_terms_estimate;
  int<lower=0> Betas_n_sigma;
  int<lower=0> Betas_n_sigma_fixed;
  int<lower=0> Betas_n_sigma_estimate;
  array[Betas_n_sigma + 1] int<lower=1, upper = Betas_raw_n_terms> Betas_re_start;
  array[Betas_n_sigma + 1] int<lower=1, upper = Betas_raw_n_terms> Betas_re_end;
  real<lower = verysmallnumber> Betas_scalarprior_sd;
  real Betas_scalarprior_mean;
  real<lower = verysmallnumber> Betas_prior_sd_sigma_estimate;
  matrix[n_geounit, Betas_raw_n_terms] Betas_model_matrix;
  matrix[Betas_raw_n_terms_fixed, k] Betas_raw_fixed;
  matrix<lower=0>[Betas_n_sigma_fixed, k] Betas_sigma_fixed;


    // // these can be used as max's when adding a level
  // // not yet used
  // real Ptilde_sigma_max;
  // real Omega_sigma_max;
  // real a_sigma_max_1;
  // real a_sigma_max_2;
  // real a_sigma_max_3;
  // real a_sigma_max_4;


    // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[fix_smoothing] real<lower=0, upper=1> Rho_fixed;
  array[fix_smoothing] real<lower=0> Tau_fixed;

