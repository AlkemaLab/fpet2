  // hierarchical parameters
  z_Omega_raw_estimate ~ std_normal();
  z_Omega_sigma_estimate ~ normal(0, z_Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(z_Epsilon_innovation) ~ std_normal();
    if(z_fix_smoothing == 0) {
      z_Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      z_Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }


