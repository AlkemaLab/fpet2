  DM1_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  DM1_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    DM1_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    DM1_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    DM1_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }
