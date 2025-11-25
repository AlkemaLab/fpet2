
  vector[z_Omega_raw_n_terms_estimate] z_Omega_raw_estimate;
  vector<lower=0>[z_Omega_n_sigma_estimate] z_Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] z_Epsilon_innovation;
  array[smoothing * (1 - z_fix_smoothing)] real<lower=0, upper=1> z_Rho_estimate;
  array[smoothing * (1 - z_fix_smoothing)] real<lower=0> z_Tau_estimate;
