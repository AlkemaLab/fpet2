  // asymptote
  vector[Ptilde_n_sigma_estimate] Ptilde_sigma_estimate = reverse(Ptilde_sigma_estimate_reverse);
  vector[Ptilde_raw_n_terms] Ptilde_star = get_mu_star(
      Ptilde_n_sigma, Ptilde_n_sigma_fixed, Ptilde_n_sigma_estimate,
      Ptilde_sigma_fixed, Ptilde_sigma_estimate,
      Ptilde_scalarprior_mean,
      Ptilde_scalarprior_sd,
      Ptilde_raw_n_terms, Ptilde_raw_n_terms_fixed, Ptilde_raw_n_terms_estimate,
      Ptilde_raw_fixed, Ptilde_raw_estimate,
      Ptilde_re_start, Ptilde_re_end);
//   vector[n_geounit] tr_Ptilde = Ptilde_low + (1 - Ptilde_low)*inv_logit(get_mu(
// use probit
   vector[n_geounit] tr_Ptilde = //rep_vector(0.9999, n_geounit);
   Ptilde_low + (1 - Ptilde_low)*Phi_approx(get_mu(
      Ptilde_star,
      Ptilde_raw_n_terms,
      n_geounit,
      Ptilde_model_matrix_w, Ptilde_model_matrix_v, Ptilde_model_matrix_u));
  // omega
  vector[Omega_raw_n_terms] Omega_star = get_mu_star(
    Omega_n_sigma, Omega_n_sigma_fixed, Omega_n_sigma_estimate,
    Omega_sigma_fixed, Omega_sigma_estimate,
    Omega_scalarprior_mean, Omega_scalarprior_sd,
    Omega_raw_n_terms, Omega_raw_n_terms_fixed, Omega_raw_n_terms_estimate,
    Omega_raw_fixed, Omega_raw_estimate,
    Omega_re_start, Omega_re_end);
  vector[n_geounit] Omega = get_mu(
    Omega_star,
    Omega_raw_n_terms,
    n_geounit,
    Omega_model_matrix_w, Omega_model_matrix_v, Omega_model_matrix_u);

  // splines coeff
  matrix<lower=0>[Betas_n_sigma_estimate, k] Betas_sigma_estimate;
  if (Betas_n_sigma_estimate > 0) {
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 1] = reverse(Betas_sigma_estimate_reverse_1);
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 2] = reverse(Betas_sigma_estimate_reverse_2);
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 3] = reverse(Betas_sigma_estimate_reverse_3);
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 4] = reverse(Betas_sigma_estimate_reverse_4);
    //Betas_sigma_estimate[1:Betas_n_sigma_estimate, 5] = reverse(Betas_sigma_estimate_reverse_5);
    //   Betas_sigma_estimate[1:Betas_n_sigma_estimate, 6] = reverse(Betas_sigma_estimate_reverse_6);
    //   Betas_sigma_estimate[1:Betas_n_sigma_estimate, 7] = reverse(Betas_sigma_estimate_reverse_7);
  }
  matrix[Betas_raw_n_terms,Betas_k_terms] Betas_star = get_mudimhk_star(k ,
       Betas_n_sigma, Betas_n_sigma_fixed, Betas_n_sigma_estimate,
       Betas_sigma_fixed, Betas_sigma_estimate,
       Betas_scalarprior_mean, Betas_scalarprior_sd,
       Betas_raw_n_terms, Betas_raw_n_terms_fixed, Betas_raw_n_terms_estimate,
       Betas_raw_fixed, Betas_raw_estimate,
       Betas_re_start, Betas_re_end);
   matrix[n_geounit,Betas_k_terms] Betas = get_mudimhk(k,
       Betas_star,
       Betas_raw_n_terms,
       n_geounit,
       Betas_model_matrix_w, Betas_model_matrix_v, Betas_model_matrix_u);
  matrix[n_geounit,Betas_k_terms] tr_Betas_nonzero = Betas_lower_bound + (Betas_upper_bound - Betas_lower_bound) * inv_logit(Betas);


 // smoothing
  matrix[n_geounit, T] Epsilon;
  Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    Rho_fixed, Rho_estimate,
    Tau_fixed, Tau_estimate,
    Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_Eta_obs = process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing, Epsilon,
                        Omega, tr_Ptilde, tr_Betas_nonzero,
                        k, num_basis,   ext_knots,  spline_degree, add_shock);
