


// RW1 model for all populations considered
matrix rw1process_model_returns_etatr(
   int n_geounit, int T, int t_star,
    array[] int t_min, array[] int t_max,
    int smoothing, matrix epsilon,
    vector Omega,
    int add_shock
   ){

  // shocks
  matrix[n_geounit, T] shock = rep_matrix(0, n_geounit, T);
  matrix[n_geounit, T] tr_eta;
  for(c in 1:n_geounit) {
    tr_eta[c, 1:T] = get_rw1_1pop_returns_etatr(Omega[c], smoothing, epsilon[c,],
            add_shock, shock[c,],
            t_star, t_min[c], t_max[c], T);
  }
  return(tr_eta);
}


// rw1 for 1 pop
row_vector get_rw1_1pop_returns_etatr(
     real Omega,
        int smoothing,
    row_vector epsilon,
    real add_shock, row_vector shock,
      int t_star, int t_min, int t_max, int T){
  row_vector[T] tr_eta;
  tr_eta[t_star] =  Omega;
   // to do: add requirement somewhere that tmax >= tstar+1
  for(t in (t_star + 1):t_max) {
    tr_eta[t] = tr_eta[t - 1]  + smoothing*epsilon[t] + add_shock*shock[t];
  }
  // to do: add requirement somewhere that tmin <= tstar-1
  for(q in 1:(t_star - t_min)) {
    int t = t_star - q;
    tr_eta[t] = tr_eta[t + 1]  - smoothing*epsilon[ t + 1] - add_shock*shock[ t+1];
  }
  // to do: fix
  // right now, added nonsense value for eta outside observation period
  // else stan crashes (can't apply the function to NAs)
  if (t_min > 1)
    tr_eta[1:(t_min-1)] = rep_row_vector(0,(t_min-1)) ;
  if (t_max < T)
    tr_eta[(t_max+1):T] = rep_row_vector(0,T - t_max) ;
  return(tr_eta);
}


// rw1 model for all populations considered outside obs period
matrix rw1process_model_outsideobs(
    matrix tr_eta, // transformed eta of dime geo_unit x T with estimates given inside obs period
    int n_geounit, int T,
    array[] int t_min, array[] int t_max,
    int smoothing, matrix epsilon,
         int add_shock
   ){
  matrix[n_geounit, T] shock = rep_matrix(0, n_geounit, T);
  matrix[n_geounit, T] tr_eta_again = tr_eta;
  for(c in 1:n_geounit) {
    tr_eta_again[c, 1:T] = get_rw1_1pop_outsideobs(
            tr_eta[c, 1:T],
            smoothing,  epsilon[c,],
            add_shock, shock[c,],
            t_min[c], t_max[c], T);
  }
  return(tr_eta_again);
}

// transition model for 1 population outside obs period
row_vector get_rw1_1pop_outsideobs(
  row_vector tr_eta_obs, // transformed eta of length T with estimates given inside obs period
     int smoothing,
    row_vector epsilon,
    // matrix epsilon,
    real add_shock, row_vector shock,
    int t_min, int t_max, int T){
  row_vector[T] tr_eta;
  tr_eta[t_min:t_max] = tr_eta_obs[t_min:t_max];
  if (t_max < T){
    for(t in (t_max + 1):T) {
      tr_eta[t] = tr_eta[t - 1]  + smoothing*epsilon[t] + add_shock*shock[t];
    }
  }
  if (t_min > 1){
    for(q in 1:(t_min - 1)) {
      int t = t_min - q;
      tr_eta[t] = tr_eta[t + 1]  - smoothing*epsilon[ t + 1] - add_shock*shock[ t+1];
    }
  }
  return(tr_eta);
}
