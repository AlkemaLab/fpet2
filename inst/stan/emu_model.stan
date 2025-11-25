  // emu data model
  if (N_emu > 0) {
   for (c in 1:n_geounit){
    log_sigma_delta[c,] ~ normal(mean_log_sigma_type, hierarchical_sigma); // one term per pop-type_index combination
   }
   for (j in 1:N_emu){
      real roc_emu_totalsd_j = sqrt(sd2_emu_roc[j] + exp(log_sigma_delta[c_emu_j[j], ss_type_index_j[j]])^2);
      // married women estimates
      real demand_jt = inv_tr_eta(tr_d_Eta_obs[c_emu_j[j], t_emu_j[j]]);
      real modern_jt = inv_tr_eta(tr_Eta_obs[c_emu_j[j], t_emu_j[j]]) * demand_jt;
      real demand_jtmin1 = inv_tr_eta(tr_d_Eta_obs[c_emu_j[j], t_emu_j[j]-1]);
      real modern_jtmin1 = inv_tr_eta(tr_Eta_obs[c_emu_j[j], t_emu_j[j]-1]) * demand_jtmin1;
      real mean_delta_emu_j;
      if (emu_for_allwomen_j[j] == 0){ // married women
        mean_delta_emu_j = modern_jt - modern_jtmin1;
      } else {
        real prop_married_jt = prop_married_rt[c_emu_j[j], t_emu_j[j]];
        real prop_married_jtmin1 = prop_married_rt[c_emu_j[j], t_emu_j[j]-1];
        // c and time are the same across married and unmarried
        real  unmarried_demand_jt = inv_tr_eta(tr_unmarried_d_Eta_obs[c_emu_j[j], t_emu_j[j]]);
        real unmarried_modern_jt = inv_tr_eta(tr_unmarried_Eta_obs[c_emu_j[j], t_emu_j[j]]) * unmarried_demand_jt;
        real allwomen_modern_jt = prop_married_jt * modern_jt +
                                              (1-prop_married_jt) * unmarried_modern_jt;
        // for tmin1
        real unmarried_demand_jtmin1 = inv_tr_eta(tr_unmarried_d_Eta_obs[c_emu_j[j], t_emu_j[j]-1]);
        real  unmarried_modern_jtmin1 = inv_tr_eta(tr_unmarried_Eta_obs[c_emu_j[j], t_emu_j[j]-1]) * unmarried_demand_jtmin1;
        real allwomen_modern_jtmin1 = prop_married_jtmin1 * modern_jtmin1 +
                                              (1-prop_married_jtmin1) * unmarried_modern_jtmin1;
        mean_delta_emu_j = allwomen_modern_jt - allwomen_modern_jtmin1;
      }
     emu_roc[j] ~ normal(mean_delta_emu_j, roc_emu_totalsd_j);
   }
  }


