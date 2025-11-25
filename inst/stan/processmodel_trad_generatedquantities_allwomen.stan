
// define trad
// earlier blocks define Eta (ds) and d_Eta (demand)
// unmet = demand - modern = (demand - ds*demand) = demand*(1-ds)
//. = d_eta * (1-Eta)
// we obtain tr_z_eta_obs, and trad = z_eta .* unmet;
// to do: we can now remove extra processing in R as unmet is already done
matrix[n_geounit, T] unmet = d_Eta.*(1-Eta);
matrix[n_geounit, T] modern = d_Eta - unmet;
matrix[n_geounit, T] trad = inv_tr_eta_matrix(
        rw1process_model_outsideobs(tr_z_Eta_obs,
                  n_geounit, T,  t_min, t_max,
                  smoothing, z_Epsilon, add_shock)).*unmet;


matrix[n_geounit, T] unmarried_unmet = unmarried_d_Eta.*(1-unmarried_Eta);
matrix[n_geounit, T] unmarried_modern = unmarried_d_Eta - unmarried_unmet;
matrix[n_geounit, T] unmarried_trad = inv_tr_eta_matrix(
        rw1process_model_outsideobs(tr_unmarried_z_Eta_obs,
                  n_geounit, T,  t_min, t_max,
                  smoothing, unmarried_z_Epsilon, add_shock)).*unmarried_unmet;
