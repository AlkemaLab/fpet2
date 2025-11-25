  // hierarchical parameters
  Ptilde_raw_estimate ~ std_normal();
  Omega_raw_estimate ~ std_normal();
  to_vector(Betas_raw_estimate) ~ std_normal();

  // variances
  //Ptilde_sigma_estimate ~ normal(0, Ptilde_prior_sd_sigma_estimate)T[0, positive_infinity()];
  Ptilde_sigma_estimate_reverse ~ normal(0, Ptilde_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  Omega_sigma_estimate ~ normal(0, Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(Betas_sigma_estimate) ~ normal(0, Betas_prior_sd_sigma_estimate)T[0, positive_infinity()];
   // one vector per spline coefficient
  to_vector(Betas_sigma_estimate_reverse_1) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(Betas_sigma_estimate_reverse_2) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(Betas_sigma_estimate_reverse_3) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(Betas_sigma_estimate_reverse_4) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(Betas_sigma_estimate_reverse_5) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(Betas_sigma_estimate_reverse_6) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(Betas_sigma_estimate_reverse_7) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(Epsilon_innovation) ~ std_normal();
    if(fix_smoothing == 0) {
      Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }


