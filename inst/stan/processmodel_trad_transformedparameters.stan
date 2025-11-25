

  // z_Omega
  vector[z_Omega_raw_n_terms] z_Omega_star = get_mu_star(
    z_Omega_n_sigma, z_Omega_n_sigma_fixed, z_Omega_n_sigma_estimate,
    z_Omega_sigma_fixed, z_Omega_sigma_estimate,
    z_Omega_scalarprior_mean, z_Omega_scalarprior_sd,
    z_Omega_raw_n_terms, z_Omega_raw_n_terms_fixed, z_Omega_raw_n_terms_estimate,
    z_Omega_raw_fixed, z_Omega_raw_estimate,
    z_Omega_re_start, z_Omega_re_end);
  vector[n_geounit] z_Omega = get_mu(
    z_Omega_star,
    z_Omega_raw_n_terms,
    n_geounit,
    z_Omega_model_matrix_w, z_Omega_model_matrix_v, z_Omega_model_matrix_u);


 // smoothing
  matrix[n_geounit, T] z_Epsilon;
  z_Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, z_fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    z_Rho_fixed, z_Rho_estimate,
    z_Tau_fixed, z_Tau_estimate,
    z_Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_z_Eta_obs = rw1process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing,
                       z_Epsilon,
                        z_Omega, add_shock);
