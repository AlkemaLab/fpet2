  matrix[n_geounit, T] Eta = inv_tr_eta_matrix(
            process_model_outsideobs(tr_Eta_obs,
                      n_geounit, T,  t_min, t_max,
                      smoothing, Epsilon,
                      tr_Ptilde, tr_Betas_nonzero,
                      k, num_basis,   ext_knots,  spline_degree, add_shock));
