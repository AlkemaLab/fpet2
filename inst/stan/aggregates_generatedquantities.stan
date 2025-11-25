  // subnational aggregates
  // assume we only do this with traditional, then we have trad and unmet
  // (and unmet, so could update this)
  // vector[T] demand_aggr;
  vector[T] modern_aggr;
  vector[T] unmet_aggr;
  vector[T] trad_aggr;
  for(t in 1:T) {
    // weighting is for the proportions (so demand, need to get modern first at regional level)
     // demand_aggr[t] = sum(d_Eta[,t] .* to_vector(geo_unit_natpop_weight_tr[t, ]));
     //vector[n_geounit] modern_regions = d_Eta[,t] .* Eta[,t];
     modern_aggr[t] = sum(modern[,t] .* to_vector(geo_unit_natpop_weight_tr[t, ]));
     unmet_aggr[t] = sum(unmet[, t] .* to_vector(geo_unit_natpop_weight_tr[t, ]));
     trad_aggr[t] = sum(trad[, t] .* to_vector(geo_unit_natpop_weight_tr[t, ]));
  }


