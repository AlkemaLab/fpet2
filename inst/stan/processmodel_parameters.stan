  // splines coefficients
  matrix[Betas_raw_n_terms_estimate, k] Betas_raw_estimate;
  // matrix<lower = 0>[a_n_sigma_estimate, k>3 ? 3 : k] a_sigma_estimate;
  // with ordered variances
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_1;
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_2;
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_3;
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_4;
  //positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_5;
  // positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_6;
  // positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_7;

  vector[Ptilde_raw_n_terms_estimate] Ptilde_raw_estimate;
  //vector<lower=verysmallnumber>[Ptilde_n_sigma_estimate] Ptilde_sigma_estimate;
  positive_ordered [Ptilde_n_sigma_estimate] Ptilde_sigma_estimate_reverse;

  vector[Omega_raw_n_terms_estimate] Omega_raw_estimate;
  vector<lower=0>[Omega_n_sigma_estimate] Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] Epsilon_innovation;
  array[smoothing * (1 - fix_smoothing)] real<lower=0, upper=1> Rho_estimate;
  array[smoothing * (1 - fix_smoothing)] real<lower=0> Tau_estimate;
