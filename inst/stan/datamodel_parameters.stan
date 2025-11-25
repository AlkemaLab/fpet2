  array[DM1_fix_nonse ? 0 : S] real<lower=0> DM1_nonse_estimate;
  array[DM1_fix_nonse ? 0 : 1] real<lower=0> DM1_sdbias_estimate;
  array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*N] DM1_local_shrinkage_dm;
  array[(1-DM1_fix_nonse)] real<lower = 0, upper = 0.99> DM1_rho_pma_estimate;
 // array[(1-DM1_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;

