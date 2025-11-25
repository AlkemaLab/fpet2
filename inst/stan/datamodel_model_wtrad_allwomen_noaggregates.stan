 // model
  // fit to the data

  // get means for all, to use for pma and non-pma
  vector[N] modern;
  vector[N] logit_unmetovernonmodern;
  vector[N] trad;
  vector[N] logit_dm3;

  for (i in 1:N){
    real demand = inv_tr_eta(tr_d_Eta_obs[geo_unit[i], time[i]]);
    modern[i] = inv_tr_eta(tr_Eta_obs[geo_unit[i], time[i]]) * demand;
    real unmet = demand - modern[i];
    logit_unmetovernonmodern[i] = logit(unmet/(1-modern[i]));

    // define trad
    // we obtain tr_z_eta_obs, and trad = z_eta .* unmet;
    trad[i] = inv_tr_eta(tr_z_Eta_obs[geo_unit[i], time[i]])*unmet;

    // for dm3, mean = ifelse(is.na(dm2), trad/unmet, trad/1minmodern)
    if (DM2_obs_isna[i] == 0){
      logit_dm3[i] = logit(trad[i]/unmet);
    } else {
      logit_dm3[i] = logit(trad[i]/(1-modern[i]));
    }
  }

  // non-pma data
  for(i in 1:N) {
    if(held_out[i] == 0 && ispma[i] == 0) {
      if(DM1_obs_isna[i] == 0) {
          DM1_y[i] ~ normal(logit(modern[i]), DM1_scale[i]);
      }
      if(DM2_obs_isna[i] == 0) {
         DM2_y[i] ~ normal(logit_unmetovernonmodern[i], DM2_scale[i]);
      }
      if(DM3_obs_isna[i] == 0) {
         DM3_y[i] ~ normal(logit_dm3[i], DM3_scale[i]);
      }
    }
  }
  // pma
  // no densities, constrained between 0 and 1
  // not yet for held out (held out just filtered out)
  // not yet for agregates
  // same for unmet
  for (index_country in pma_country_indices){
    array[npma_c[index_country]] int indices_inc = indices_pma[nstart_pma_c[index_country]:
          (nstart_pma_c[index_country]+npma_c[index_country] - 1)];
    matrix[npma_c[index_country], npma_c[index_country]] DM1_cov_matrix_pma;
    matrix[npma_c[index_country], npma_c[index_country]] DM2_cov_matrix_pma;
    matrix[npma_c[index_country], npma_c[index_country]] DM3_cov_matrix_pma;
    for (j1 in 1:npma_c[index_country]){
      for (j2 in 1:npma_c[index_country]){
        int j1_index = indices_inc[j1];
        int j2_index = indices_inc[j2];
        DM1_cov_matrix_pma[j1, j2] =
            DM1_scale[j1_index]* DM1_scale[j2_index] *
              DM1_rho_pma^(abs(exact_time[j1_index] - exact_time[j2_index]));
        DM2_cov_matrix_pma[j1, j2] =
            DM2_scale[j1_index]* DM2_scale[j2_index] *
              DM2_rho_pma^(abs(exact_time[j1_index] - exact_time[j2_index]));
        DM3_cov_matrix_pma[j1, j2] =
          DM3_scale[j1_index]* DM3_scale[j2_index] *
            DM3_rho_pma^(abs(exact_time[j1_index] - exact_time[j2_index]));

    }}
    DM1_y[indices_inc] ~ multi_normal(logit(modern[indices_inc]), DM1_cov_matrix_pma);
    DM2_y[indices_inc] ~ multi_normal(logit_unmetovernonmodern[indices_inc], DM2_cov_matrix_pma);
    DM3_y[indices_inc] ~ multi_normal(logit_dm3[indices_inc], DM3_cov_matrix_pma);
  }



  // get means for all, to use for pma and non-pma
  vector[unmarried_N] unmarried_modern;
  vector[unmarried_N] unmarried_logit_unmetovernonmodern;
  vector[unmarried_N] unmarried_trad;
  vector[unmarried_N] unmarried_logit_dm3;

  for (i in 1:unmarried_N){
    real unmarried_demand = inv_tr_eta(tr_unmarried_d_Eta_obs[unmarried_geo_unit[i], unmarried_time[i]]);
    unmarried_modern[i] = inv_tr_eta(tr_unmarried_Eta_obs[unmarried_geo_unit[i], unmarried_time[i]]) *
                                        unmarried_demand;
    real unmarried_unmet = unmarried_demand - unmarried_modern[i];
    unmarried_logit_unmetovernonmodern[i] = logit(unmarried_unmet/(1-unmarried_modern[i]));

    // define trad
    // we obtain tr_z_eta_obs, and trad = z_eta .* unmet;
    unmarried_trad[i] = inv_tr_eta(tr_unmarried_z_Eta_obs[unmarried_geo_unit[i], unmarried_time[i]])*unmarried_unmet;

    // for dm3, mean = ifelse(is.na(dm2), trad/unmet, trad/1minmodern)
    if (unmarried_DM2_obs_isna[i] == 0){
      unmarried_logit_dm3[i] = logit(unmarried_trad[i]/unmarried_unmet);
    } else {
      unmarried_logit_dm3[i] = logit(unmarried_trad[i]/(1-unmarried_modern[i]));
    }
  }

  // non-pma data
  for(i in 1:unmarried_N) {
    if(unmarried_held_out[i] == 0 && unmarried_ispma[i] == 0) {
      if(unmarried_DM1_obs_isna[i] == 0) {
          unmarried_DM1_y[i] ~ normal(logit(unmarried_modern[i]), unmarried_DM1_scale[i]);
      }
      if(unmarried_DM2_obs_isna[i] == 0) {
         unmarried_DM2_y[i] ~ normal(unmarried_logit_unmetovernonmodern[i], unmarried_DM2_scale[i]);
      }
      if(unmarried_DM3_obs_isna[i] == 0) {
         unmarried_DM3_y[i] ~ normal(unmarried_logit_dm3[i], unmarried_DM3_scale[i]);
      }
    }
  }
  // pma
  // no densities, constrained between 0 and 1
  // not yet for held out (held out just filtered out)
  // not yet for agregates
  // same for unmet
  for (index_country in unmarried_pma_country_indices){
    array[unmarried_npma_c[index_country]] int unmarried_indices_inc =
                unmarried_indices_pma[unmarried_nstart_pma_c[index_country]:
          (unmarried_nstart_pma_c[index_country] + unmarried_npma_c[index_country] - 1)];
    matrix[unmarried_npma_c[index_country], unmarried_npma_c[index_country]] unmarried_DM1_cov_matrix_pma;
    matrix[unmarried_npma_c[index_country], unmarried_npma_c[index_country]] unmarried_DM2_cov_matrix_pma;
    matrix[unmarried_npma_c[index_country], unmarried_npma_c[index_country]] unmarried_DM3_cov_matrix_pma;
    for (j1 in 1:unmarried_npma_c[index_country]){
      for (j2 in 1:unmarried_npma_c[index_country]){
        int unmarried_j1_index = unmarried_indices_inc[j1];
        int unmarried_j2_index = unmarried_indices_inc[j2];
        unmarried_DM1_cov_matrix_pma[j1, j2] =
            unmarried_DM1_scale[unmarried_j1_index]* unmarried_DM1_scale[unmarried_j2_index] *
              unmarried_DM1_rho_pma^(abs(unmarried_exact_time[unmarried_j1_index] -
              unmarried_exact_time[unmarried_j2_index]));
        unmarried_DM2_cov_matrix_pma[j1, j2] =
            unmarried_DM2_scale[unmarried_j1_index]* unmarried_DM2_scale[unmarried_j2_index] *
              unmarried_DM2_rho_pma^(abs(unmarried_exact_time[unmarried_j1_index] -
              unmarried_exact_time[unmarried_j2_index]));
        unmarried_DM3_cov_matrix_pma[j1, j2] =
          unmarried_DM3_scale[unmarried_j1_index]* unmarried_DM3_scale[unmarried_j2_index] *
            unmarried_DM3_rho_pma^(abs(unmarried_exact_time[unmarried_j1_index] -
            unmarried_exact_time[unmarried_j2_index]));

    }}
    unmarried_DM1_y[unmarried_indices_inc] ~ multi_normal(logit(unmarried_modern[unmarried_indices_inc]),
                      unmarried_DM1_cov_matrix_pma);
    unmarried_DM2_y[unmarried_indices_inc] ~ multi_normal(unmarried_logit_unmetovernonmodern[unmarried_indices_inc],
                      unmarried_DM2_cov_matrix_pma);
    unmarried_DM3_y[unmarried_indices_inc] ~ multi_normal(unmarried_logit_dm3[unmarried_indices_inc],
                      unmarried_DM3_cov_matrix_pma);
  }
