
// transition model for all populations considered
matrix process_model_returns_etatr(
    int n_geounit, int T, int t_star,
    array[] int t_min, array[] int t_max,
    int smoothing, matrix epsilon,
    vector Omega,
    vector Ptilde,
    matrix tr_a_nonzero,
     int k, int num_basis,  vector ext_knots, int spline_degree,
         int add_shock
   ){

  matrix[n_geounit, num_basis] a;
  // set last num_constrained_zero coefficients to 0
  for(c in 1:n_geounit) {
     a[c, (k + 1):num_basis] = rep_row_vector(0, num_basis - k);
  }
  a[ , 1:k] = tr_a_nonzero;
  matrix[n_geounit, T] shock = rep_matrix(0, n_geounit, T);
  matrix[n_geounit, T] tr_eta;
  for(c in 1:n_geounit) {
    tr_eta[c, 1:T] = get_transition_1pop_returns_etatr(a[c,], Omega[c], Ptilde[c],
            smoothing,  epsilon[c,],
            add_shock, shock[c,],
            //  int t_star, int t_min, int t_max, int T,
            t_star, t_min[c], t_max[c], T,
            ext_knots, num_basis, spline_degree);
  }
  return(tr_eta);
}


// transition model for 1 population
row_vector get_transition_1pop_returns_etatr(
    row_vector a,  real Omega, real P_tilde,
     int smoothing,
    row_vector epsilon,
    // matrix epsilon,
    real add_shock, row_vector shock,
    int t_star, int t_min, int t_max, int T,
    vector ext_knots, int num_basis, int spline_degree){
  row_vector[T] tr_eta;
  real transition_function;
  tr_eta[t_star] =  Omega;
  // to do: add requirement somewhere that tmax >= tstar+1
  for(t in (t_star + 1):t_max) {
    transition_function =
      rate_spline(
       inv_tr_eta(tr_eta[t - 1]),
       P_tilde,
       a, ext_knots, num_basis, spline_degree);
    tr_eta[t] = tr_eta[t - 1] + transition_function + smoothing*epsilon[t] + add_shock*shock[t];
  }
  // to do: add requirement somewhere that tmin <= tstar-1
  for(q in 1:(t_star - t_min)) {
    int t = t_star - q;
    transition_function =
    rate_spline(
        //0;
    //     rate_spline_noasymptote(
    // //    rate_constant(
      inv_tr_eta(tr_eta[t + 1]),
      P_tilde,
      a, ext_knots, num_basis, spline_degree);
    tr_eta[t] = tr_eta[t + 1] - transition_function - smoothing*epsilon[ t + 1] - add_shock*shock[ t+1];
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


// transition model for all populations considered outside obs period
matrix process_model_outsideobs(
  matrix tr_eta, // transformed eta of dime geo_unit x T with estimates given inside obs period
    int n_geounit, int T,
    array[] int t_min, array[] int t_max,
    int smoothing, matrix epsilon,
    vector Ptilde,
    matrix tr_a_nonzero,
     int k, int num_basis,  vector ext_knots, int spline_degree,
         int add_shock
   ){

  matrix[n_geounit, num_basis] a;
  // set last num_constrained_zero coefficients to 0
  for(c in 1:n_geounit) {
     a[c, (k + 1):num_basis] = rep_row_vector(0, num_basis - k);
  }
  a[ , 1:k] = tr_a_nonzero;
  matrix[n_geounit, T] shock = rep_matrix(0, n_geounit, T);
  matrix[n_geounit, T] tr_eta_again = tr_eta;
  for(c in 1:n_geounit) {
    tr_eta_again[c, 1:T] = get_transition_1pop_outsideobs(
            tr_eta[c, 1:T],
            a[c,],  Ptilde[c],
            smoothing,  epsilon[c,],
            add_shock, shock[c,],
            t_min[c], t_max[c], T,
            ext_knots, num_basis, spline_degree);
  }
  return(tr_eta_again);
}

// transition model for 1 population outside obs period
row_vector get_transition_1pop_outsideobs(
  row_vector tr_eta_obs, // transformed eta of length T with estimates given inside obs period
    row_vector a,  real P_tilde,
     int smoothing,
    row_vector epsilon,
    // matrix epsilon,
    real add_shock, row_vector shock,
    int t_min, int t_max, int T,
    vector ext_knots, int num_basis, int spline_degree){
  real transition_function;
  row_vector[T] tr_eta;
  tr_eta[t_min:t_max] = tr_eta_obs[t_min:t_max];
  if (t_max < T){
    for(t in (t_max + 1):T) {
      transition_function = rate_spline(inv_tr_eta(tr_eta[t - 1]),
       P_tilde,
       a, ext_knots, num_basis, spline_degree);
      tr_eta[t] = tr_eta[t - 1] + transition_function + smoothing*epsilon[t] + add_shock*shock[t];
    }
  }
  if (t_min > 1){
    for(q in 1:(t_min - 1)) {
      int t = t_min - q;
      transition_function = rate_spline(inv_tr_eta(tr_eta[t + 1]),
        P_tilde,
        a, ext_knots, num_basis, spline_degree);
      tr_eta[t] = tr_eta[t + 1] - transition_function - smoothing*epsilon[ t + 1] - add_shock*shock[ t+1];
    }
  }
  return(tr_eta);
}

// transformation functions used in process model to avoid extreme outcomes
real inv_tr_eta(real x) {
  return Phi_approx(x);

}

row_vector inv_tr_eta_vector(row_vector x) {
   return Phi_approx(x);
}

vector inv_tr_eta_colvector(vector x) {
   return Phi_approx(x);
}

matrix inv_tr_eta_matrix(matrix x) {
   return Phi_approx(x);
}

row_vector tr_eta_vector(row_vector y) {
   return inv_Phi(y);
}

vector tr_eta_colvector(vector y) {
   return inv_Phi(y);
}

real rate_spline(real P, real P_tilde, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
// for testing without Ptilde (if to be used, update to remove Ptilde from arguments)
//      return deboor(P, ext_knots, a, spline_degree);
      return deboor(P / P_tilde, ext_knots, a, spline_degree);
    }

//   // extra functions for simplified fitting, non-efficient implementation to follow same for as rate_spline
//   // not yet fully finished/tested
//   real rate_constant(real P, real P_tilde, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
//     return a[1]; // * P / P_tilde;
//     }
//
//
//  real rate_spline_noasymptote(real P, real P_tilde, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
// //      row_vector[num_basis] spline_coeff = rep_row_vector(a[1], num_basis);
//       row_vector[num_basis] spline_coeff = rep_row_vector(0, num_basis);
//     //   vector[num_basis] weights = rep_vector(1, num_basis);
//     // return sum(weights*spline_coeff);
//     return deboor(P, ext_knots, spline_coeff, spline_degree);
//   }
//
//
//   real rate_constant_asymptote(real P, real P_tilde, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
//     row_vector[num_basis] spline_coeff = rep_row_vector(a[1], num_basis);
//     return deboor(P / P_tilde, ext_knots, spline_coeff, spline_degree);
//   }





