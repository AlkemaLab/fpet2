  // subnational estimation:
  // aggregation, to finish
   matrix[n_agg_units, T] tr_eta_agg;
  // if (!subnational){
  //   // not a subnat run, so not used. fix at value between 0 and 1
    for(aggunit in 1:n_agg_units) {
      for(t in 1:T) {
        tr_eta_agg[aggunit,t] = 0; // use eta_obs
      }
    }
  // } else {
  //   // TO DO: finish
  //   for(i in 1:n_agg_units) {
  //     for(t in 1:T) {
  //       eta_agg[i, t] = sum(eta[, t] .* to_vector(geo_unit_pop_wt[i, , t]));
  //     }
  //   }
  // }


