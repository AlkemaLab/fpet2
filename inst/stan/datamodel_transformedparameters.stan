  // minor note that this includes obs that are NA, could subset instead
   array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_caux_dm_estimate = DM1_sqrt_caux_dm_estimate^2;
   vector[N] DM1_scale;
   if (add_dataoutliers){
    DM1_scale = get_scale(
        DM1_local_shrinkage_dm, DM1_nonse_fixed, DM1_global_shrinkage_dm_fixed, DM1_caux_dm_fixed,
        DM1_nonse_estimate, DM1_global_shrinkage_dm_estimate, DM1_caux_dm_estimate,
        N, DM1_s, isdhs, source, DM1_nooutlier, S, DM1_fix_nonse,
        DM1_any_bias, DM1_sdbias_fixed, DM1_sdbias_estimate);
   } else {
    DM1_scale = get_scale_nooutliers(
        DM1_nonse_fixed,
        DM1_nonse_estimate,
        N, DM1_s, isdhs, source, DM1_nooutlier, S, DM1_fix_nonse);
   }
     // pma
  real<lower = 0, upper = 1> DM1_rho_pma;
  if (DM1_fix_nonse) {
    DM1_rho_pma = DM1_rho_pma_fixed[1];
  } else {
    DM1_rho_pma = DM1_rho_pma_estimate[1];
  }
  // real<lower = 0, upper = 1> z_rho_pma;
  // if (add_z && fix_z) {
  //   z_rho_pma = z_rho_pma_fixed[1];
  // } else {
  //   z_rho_pma = z_rho_pma_estimate[1];
  // }

