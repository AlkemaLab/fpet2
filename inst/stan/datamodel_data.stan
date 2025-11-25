  // one set per dm
  int<lower=0, upper=1> DM1_fix_nonse; // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  vector[N] DM1_y;                // Observations
  vector[N] DM1_s;                // Standard deviation
  array[N] int<lower=0> DM1_obs_isna;
  array[N] int<lower=0> DM1_any_bias;
  array[N] int<lower=0> DM1_nooutlier;                // should outlier be included?
  array[DM1_fix_nonse ? S : 0] real<lower=0> DM1_nonse_fixed;
  array[add_dataoutliers*DM1_fix_nonse ? 1 : 0] real<lower=0> DM1_global_shrinkage_dm_fixed;
  array[add_dataoutliers*DM1_fix_nonse ? 1 : 0] real<lower=0> DM1_caux_dm_fixed;
  array[DM1_fix_nonse ? 1 : 0] real<lower=0> DM1_sdbias_fixed;
  array[DM1_fix_nonse ? 1 : 0] real<lower=0, upper = 1> DM1_rho_pma_fixed;
  //array[DM1_fix_nonse*add_z*fix_z ? 1 : 0] real<lower=0, upper = 1> z_rho_pma_fixed;
