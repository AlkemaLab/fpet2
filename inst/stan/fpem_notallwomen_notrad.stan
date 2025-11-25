// Note on updating stan model files:
// fpem_buildingblocks.stan is the stan model used in `write_model`
// to produce stan models for the different use cases.
// So any updates need to made in fpem_buildingblocks.stan
// Then write_model can be called to get the corresponding stan models.

functions {

#include ./function_deboor.stan
#include ./functions_dm.stan
#include ./function_eps.stan
#include ./functions_hierarchical_params.stan
#include ./functions_processmodel.stan
#include ./functions_rwprocessmodel.stan
#include ./functions_truncatednormal.stan

}
/////////////////////////////////////////////////////
data {
  // same regardless of number of indicators/data models//martial group or all
   real<lower = 0> verysmallnumber; // lower bound for sigmas
  int add_shock; // doesn't work yet
  int add_dataoutliers;
  int generate_quantities;
  int validation_run;
  int T; // Number of time points
  int<lower=1, upper=T> t_star;
  // we need to have n_geounit and t_min and t_max the same across
  int n_geounit; // Number of lowest-level geographic units, e.g. countries or states
  array[n_geounit] int<lower=1, upper=T> t_min;
  array[n_geounit] int<lower=1, upper=T> t_max;

  // settings, same across indicators
  // process model
  // smoothing
  // 0: no smoothing component. 1: AR(1) smoothing component.
  int<lower=0, upper=1> smoothing;

  // for demand and demand satisfied, not for traditional
  int<lower=0, upper=1> fix_smoothing;// 0: don't fix smoothing component. 1: fix AR(1) smoothing component
  // if fixed, we fix rho and tau.
  // data model
  // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  // int<lower=0, upper=1> DM2_fix_nonse;
  // int<lower=0, upper=1> DM2_fix_nonse;
  // if fix_nonse == 1, one nonse standard deviation estimate per source type

  // related to splines (some for previous use in generated quantities)
  // currently we keep these the same across marital groups
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  // int num_constrained_zero;
  int num_basis;
  int k; // number of non-zero splines coefficients
  int num_grid;
  vector[num_grid] grid;
  matrix[num_knots + spline_degree - 1, num_grid] B;

    // data properties, different between marital groups
  int N; // Number of observations
  int S; // Number of sources
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=0, upper=n_geounit> geo_unit; // Geographic unit of each observation, 0 if observation is for an aggregated unit
  // updated to lower of 0 for subnational runs
  array[N] int<lower=0, upper=S> source;   // Source of each observation
  array[N] int<lower=0, upper=1> held_out; // Whether to hold out each observation
  array[N] int<lower=0> isdhs;            // DHS?
  // pma data
  array[N] real<lower=0> ispma;            // PMA?
  array[N] real<lower=0> exact_time;
  int<lower=0> N_pma_country_indices;
  array[N_pma_country_indices] int<lower=1, upper = n_geounit> pma_country_indices;
  array[n_geounit] int npma_c;
  array[n_geounit] int nstart_pma_c;
  int<lower=0> N_pma_obs_indices;
  array[N_pma_obs_indices] int indices_pma;


  // for parameters of process model
  // in old code, "d_" referred to demand
  // here do ds
  //modern[m] = eta[m] .* d_eta[m];
  //unmet[m] = d_eta[m] - modern[m];

    real Betas_lower_bound;
  real Betas_upper_bound;
  real Ptilde_low;

  int<lower=0> Ptilde_raw_n_terms;
  int<lower=0> Ptilde_raw_n_terms_fixed;
  int<lower=0> Ptilde_raw_n_terms_estimate;
  int<lower=0> Ptilde_n_sigma;
  int<lower=0> Ptilde_n_sigma_fixed;
  int<lower=0> Ptilde_n_sigma_estimate;
  array[Ptilde_n_sigma + 1] int<lower=1, upper=Ptilde_raw_n_terms> Ptilde_re_start;
  array[Ptilde_n_sigma + 1] int<lower=1, upper=Ptilde_raw_n_terms> Ptilde_re_end;
  matrix[n_geounit, Ptilde_raw_n_terms] Ptilde_model_matrix;
  real<lower = verysmallnumber> Ptilde_scalarprior_sd;
  real Ptilde_scalarprior_mean;
  real<lower = verysmallnumber> Ptilde_prior_sd_sigma_estimate;
  vector[Ptilde_raw_n_terms_fixed] Ptilde_raw_fixed;
  vector<lower=0>[Ptilde_n_sigma_fixed] Ptilde_sigma_fixed;

  int<lower=0> Omega_raw_n_terms;
  int<lower=0> Omega_n_sigma;
  int<lower=0> Omega_raw_n_terms_fixed;
  int<lower=0> Omega_raw_n_terms_estimate;
  int<lower=0> Omega_n_sigma_fixed;
  int<lower=0> Omega_n_sigma_estimate;
  array[Omega_n_sigma + 1] int<lower=1, upper=Omega_raw_n_terms> Omega_re_start;
  array[Omega_n_sigma + 1] int<lower=1, upper=Omega_raw_n_terms> Omega_re_end;
  matrix[n_geounit, Omega_raw_n_terms] Omega_model_matrix;
  real<lower = verysmallnumber> Omega_scalarprior_sd;
  real Omega_scalarprior_mean;
  real<lower = verysmallnumber> Omega_prior_sd_sigma_estimate;
  vector[Omega_raw_n_terms_fixed] Omega_raw_fixed;
  vector<lower=0>[Omega_n_sigma_fixed] Omega_sigma_fixed;

  int Betas_k_terms; // same as k; to do clean to have just 1
  int<lower=0> Betas_raw_n_terms;
  int<lower=0> Betas_raw_n_terms_fixed;
  int<lower=0> Betas_raw_n_terms_estimate;
  int<lower=0> Betas_n_sigma;
  int<lower=0> Betas_n_sigma_fixed;
  int<lower=0> Betas_n_sigma_estimate;
  array[Betas_n_sigma + 1] int<lower=1, upper = Betas_raw_n_terms> Betas_re_start;
  array[Betas_n_sigma + 1] int<lower=1, upper = Betas_raw_n_terms> Betas_re_end;
  real<lower = verysmallnumber> Betas_scalarprior_sd;
  real Betas_scalarprior_mean;
  real<lower = verysmallnumber> Betas_prior_sd_sigma_estimate;
  matrix[n_geounit, Betas_raw_n_terms] Betas_model_matrix;
  matrix[Betas_raw_n_terms_fixed, k] Betas_raw_fixed;
  matrix<lower=0>[Betas_n_sigma_fixed, k] Betas_sigma_fixed;


    // // these can be used as max's when adding a level
  // // not yet used
  // real Ptilde_sigma_max;
  // real Omega_sigma_max;
  // real a_sigma_max_1;
  // real a_sigma_max_2;
  // real a_sigma_max_3;
  // real a_sigma_max_4;


    // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[fix_smoothing] real<lower=0, upper=1> Rho_fixed;
  array[fix_smoothing] real<lower=0> Tau_fixed;


    real d_Betas_lower_bound;
  real d_Betas_upper_bound;
  real d_Ptilde_low;

  int<lower=0> d_Ptilde_raw_n_terms;
  int<lower=0> d_Ptilde_raw_n_terms_fixed;
  int<lower=0> d_Ptilde_raw_n_terms_estimate;
  int<lower=0> d_Ptilde_n_sigma;
  int<lower=0> d_Ptilde_n_sigma_fixed;
  int<lower=0> d_Ptilde_n_sigma_estimate;
  array[d_Ptilde_n_sigma + 1] int<lower=1, upper=d_Ptilde_raw_n_terms> d_Ptilde_re_start;
  array[d_Ptilde_n_sigma + 1] int<lower=1, upper=d_Ptilde_raw_n_terms> d_Ptilde_re_end;
  matrix[n_geounit, d_Ptilde_raw_n_terms] d_Ptilde_model_matrix;
  real<lower = verysmallnumber> d_Ptilde_scalarprior_sd;
  real d_Ptilde_scalarprior_mean;
  real<lower = verysmallnumber> d_Ptilde_prior_sd_sigma_estimate;
  vector[d_Ptilde_raw_n_terms_fixed] d_Ptilde_raw_fixed;
  vector<lower=0>[d_Ptilde_n_sigma_fixed] d_Ptilde_sigma_fixed;

  int<lower=0> d_Omega_raw_n_terms;
  int<lower=0> d_Omega_n_sigma;
  int<lower=0> d_Omega_raw_n_terms_fixed;
  int<lower=0> d_Omega_raw_n_terms_estimate;
  int<lower=0> d_Omega_n_sigma_fixed;
  int<lower=0> d_Omega_n_sigma_estimate;
  array[d_Omega_n_sigma + 1] int<lower=1, upper=d_Omega_raw_n_terms> d_Omega_re_start;
  array[d_Omega_n_sigma + 1] int<lower=1, upper=d_Omega_raw_n_terms> d_Omega_re_end;
  matrix[n_geounit, d_Omega_raw_n_terms] d_Omega_model_matrix;
  real<lower = verysmallnumber> d_Omega_scalarprior_sd;
  real d_Omega_scalarprior_mean;
  real<lower = verysmallnumber> d_Omega_prior_sd_sigma_estimate;
  vector[d_Omega_raw_n_terms_fixed] d_Omega_raw_fixed;
  vector<lower=0>[d_Omega_n_sigma_fixed] d_Omega_sigma_fixed;

  int d_Betas_k_terms; // same as k; to do clean to have just 1
  int<lower=0> d_Betas_raw_n_terms;
  int<lower=0> d_Betas_raw_n_terms_fixed;
  int<lower=0> d_Betas_raw_n_terms_estimate;
  int<lower=0> d_Betas_n_sigma;
  int<lower=0> d_Betas_n_sigma_fixed;
  int<lower=0> d_Betas_n_sigma_estimate;
  array[d_Betas_n_sigma + 1] int<lower=1, upper = d_Betas_raw_n_terms> d_Betas_re_start;
  array[d_Betas_n_sigma + 1] int<lower=1, upper = d_Betas_raw_n_terms> d_Betas_re_end;
  real<lower = verysmallnumber> d_Betas_scalarprior_sd;
  real d_Betas_scalarprior_mean;
  real<lower = verysmallnumber> d_Betas_prior_sd_sigma_estimate;
  matrix[n_geounit, d_Betas_raw_n_terms] d_Betas_model_matrix;
  matrix[d_Betas_raw_n_terms_fixed, k] d_Betas_raw_fixed;
  matrix<lower=0>[d_Betas_n_sigma_fixed, k] d_Betas_sigma_fixed;


    // // these can be used as max's when adding a level
  // // not yet used
  // real d_Ptilde_sigma_max;
  // real d_Omega_sigma_max;
  // real a_sigma_max_1;
  // real a_sigma_max_2;
  // real a_sigma_max_3;
  // real a_sigma_max_4;


    // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[fix_smoothing] real<lower=0, upper=1> d_Rho_fixed;
  array[fix_smoothing] real<lower=0> d_Tau_fixed;


   



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

    // one set per dm
  int<lower=0, upper=1> DM2_fix_nonse; // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  vector[N] DM2_y;                // Observations
  vector[N] DM2_s;                // Standard deviation
  array[N] int<lower=0> DM2_obs_isna;
  array[N] int<lower=0> DM2_any_bias;
  array[N] int<lower=0> DM2_nooutlier;                // should outlier be included?
  array[DM2_fix_nonse ? S : 0] real<lower=0> DM2_nonse_fixed;
  array[add_dataoutliers*DM2_fix_nonse ? 1 : 0] real<lower=0> DM2_global_shrinkage_dm_fixed;
  array[add_dataoutliers*DM2_fix_nonse ? 1 : 0] real<lower=0> DM2_caux_dm_fixed;
  array[DM2_fix_nonse ? 1 : 0] real<lower=0> DM2_sdbias_fixed;
  array[DM2_fix_nonse ? 1 : 0] real<lower=0, upper = 1> DM2_rho_pma_fixed;
  //array[DM2_fix_nonse*add_z*fix_z ? 1 : 0] real<lower=0, upper = 1> z_rho_pma_fixed;

   

   

  // doesn't do anything
  // subnational
  int n_agg_units;                            // Number of geo units observed at aggregated spatial scale (e.g. number of nations)
  array[N] int<lower=0, upper=1> is_agg_obs; // Whether each observation is for a lowest-level geo unit or aggregated geo units
  array[N] int<lower=0, upper=n_agg_units> agg_unit; // For observations at aggregated spatial scales, index of aggregated unit, 0 if observation is not for an aggregated unit
  array[n_agg_units, n_geounit, T] real geo_unit_pop_wt; // Weights used to aggregate geo units. Each column should sum to 1
  // subnational indicator determines whether aggregates are calculated or not
  int<lower=0, upper=1> subnational;
  // 0: smoothing terms independent across lowest-level geo units
  // 1: smoothing terms correlated among lowest-level geo units in the same group
  int<lower=0, upper=1> correlated_smoothing;
  int n_cor_smoothing_blocks; // basically, number of national level units
  // number of low-level geo units in each block
  array[n_cor_smoothing_blocks] int<lower=1, upper=n_geounit> cor_smoothing_block_sizes;
  int max_cor_smoothing_block_size;

  int<lower=0, upper=1> fix_subnat_corr;
  array[fix_subnat_corr ? 1 : 0]  real<lower=0, upper = 1> rho_correlationeps_fixed;

      

}

/////////////////////////////////////////////////////
transformed data {
  // splines set up: ext_knots (could pass from R too)
  //int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  //int num_constrained_zero = spline_degree + 1;
  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;

  // transformations for hier parameters
   vector[rows(csr_extract_w(Ptilde_model_matrix))] Ptilde_model_matrix_w    = csr_extract_w(Ptilde_model_matrix);
  array[size(csr_extract_v(Ptilde_model_matrix))] int Ptilde_model_matrix_v = csr_extract_v(Ptilde_model_matrix);
  array[size(csr_extract_u(Ptilde_model_matrix))] int Ptilde_model_matrix_u = csr_extract_u(Ptilde_model_matrix);

  vector[rows(csr_extract_w(Omega_model_matrix))] Omega_model_matrix_w     = csr_extract_w(Omega_model_matrix);
  array[size(csr_extract_v(Omega_model_matrix))] int Omega_model_matrix_v  = csr_extract_v(Omega_model_matrix);
  array[size(csr_extract_u(Omega_model_matrix))] int Omega_model_matrix_u  = csr_extract_u(Omega_model_matrix);

  vector[rows(csr_extract_w(Betas_model_matrix))] Betas_model_matrix_w    = csr_extract_w(Betas_model_matrix);
  array[size(csr_extract_v(Betas_model_matrix))] int Betas_model_matrix_v = csr_extract_v(Betas_model_matrix);
  array[size(csr_extract_u(Betas_model_matrix))] int Betas_model_matrix_u = csr_extract_u(Betas_model_matrix);

   vector[rows(csr_extract_w(d_Ptilde_model_matrix))] d_Ptilde_model_matrix_w    = csr_extract_w(d_Ptilde_model_matrix);
  array[size(csr_extract_v(d_Ptilde_model_matrix))] int d_Ptilde_model_matrix_v = csr_extract_v(d_Ptilde_model_matrix);
  array[size(csr_extract_u(d_Ptilde_model_matrix))] int d_Ptilde_model_matrix_u = csr_extract_u(d_Ptilde_model_matrix);

  vector[rows(csr_extract_w(d_Omega_model_matrix))] d_Omega_model_matrix_w     = csr_extract_w(d_Omega_model_matrix);
  array[size(csr_extract_v(d_Omega_model_matrix))] int d_Omega_model_matrix_v  = csr_extract_v(d_Omega_model_matrix);
  array[size(csr_extract_u(d_Omega_model_matrix))] int d_Omega_model_matrix_u  = csr_extract_u(d_Omega_model_matrix);

  vector[rows(csr_extract_w(d_Betas_model_matrix))] d_Betas_model_matrix_w    = csr_extract_w(d_Betas_model_matrix);
  array[size(csr_extract_v(d_Betas_model_matrix))] int d_Betas_model_matrix_v = csr_extract_v(d_Betas_model_matrix);
  array[size(csr_extract_u(d_Betas_model_matrix))] int d_Betas_model_matrix_u = csr_extract_u(d_Betas_model_matrix);

   


}

/////////////////////////////////////////////////////
parameters {

  // for parameters of process model
    // splines coefficients
  matrix[Betas_raw_n_terms_estimate, k] Betas_raw_estimate;
  // matrix<lower = 0>[a_n_sigma_estimate, k>3 ? 3 : k] a_sigma_estimate;
  // with ordered variances
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_1;
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_2;
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_3;
  positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_4;
  //positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_5;
  // positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_6;
  // positive_ordered [Betas_n_sigma_estimate] Betas_sigma_estimate_reverse_7;

  vector[Ptilde_raw_n_terms_estimate] Ptilde_raw_estimate;
  //vector<lower=verysmallnumber>[Ptilde_n_sigma_estimate] Ptilde_sigma_estimate;
  positive_ordered [Ptilde_n_sigma_estimate] Ptilde_sigma_estimate_reverse;

  vector[Omega_raw_n_terms_estimate] Omega_raw_estimate;
  vector<lower=0>[Omega_n_sigma_estimate] Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] Epsilon_innovation;
  array[smoothing * (1 - fix_smoothing)] real<lower=0, upper=1> Rho_estimate;
  array[smoothing * (1 - fix_smoothing)] real<lower=0> Tau_estimate;

    // splines coefficients
  matrix[d_Betas_raw_n_terms_estimate, k] d_Betas_raw_estimate;
  // matrix<lower = 0>[a_n_sigma_estimate, k>3 ? 3 : k] a_sigma_estimate;
  // with ordered variances
  positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_1;
  positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_2;
  positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_3;
  positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_4;
  //positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_5;
  // positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_6;
  // positive_ordered [d_Betas_n_sigma_estimate] d_Betas_sigma_estimate_reverse_7;

  vector[d_Ptilde_raw_n_terms_estimate] d_Ptilde_raw_estimate;
  //vector<lower=verysmallnumber>[d_Ptilde_n_sigma_estimate] d_Ptilde_sigma_estimate;
  positive_ordered [d_Ptilde_n_sigma_estimate] d_Ptilde_sigma_estimate_reverse;

  vector[d_Omega_raw_n_terms_estimate] d_Omega_raw_estimate;
  vector<lower=0>[d_Omega_n_sigma_estimate] d_Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] d_Epsilon_innovation;
  array[smoothing * (1 - fix_smoothing)] real<lower=0, upper=1> d_Rho_estimate;
  array[smoothing * (1 - fix_smoothing)] real<lower=0> d_Tau_estimate;

    

  // Data model, one set per dm
    array[DM1_fix_nonse ? 0 : S] real<lower=0> DM1_nonse_estimate;
  array[DM1_fix_nonse ? 0 : 1] real<lower=0> DM1_sdbias_estimate;
  array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*N] DM1_local_shrinkage_dm;
  array[(1-DM1_fix_nonse)] real<lower = 0, upper = 0.99> DM1_rho_pma_estimate;
 // array[(1-DM1_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


    array[DM2_fix_nonse ? 0 : S] real<lower=0> DM2_nonse_estimate;
  array[DM2_fix_nonse ? 0 : 1] real<lower=0> DM2_sdbias_estimate;
  array[add_dataoutliers * (1-DM2_fix_nonse)] real<lower=0> DM2_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-DM2_fix_nonse)] real<lower=0> DM2_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*N] DM2_local_shrinkage_dm;
  array[(1-DM2_fix_nonse)] real<lower = 0, upper = 0.99> DM2_rho_pma_estimate;
 // array[(1-DM2_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


    

  //{{subnational_parameters}}
  // for subnational
  array[fix_subnat_corr ? 0 : 1] real<lower=0, upper = 1> rho_correlationeps_estimate;   // for correlated eps
    

}

/////////////////////////////////////////////////////
transformed parameters {

    // asymptote
  vector[Ptilde_n_sigma_estimate] Ptilde_sigma_estimate = reverse(Ptilde_sigma_estimate_reverse);
  vector[Ptilde_raw_n_terms] Ptilde_star = get_mu_star(
      Ptilde_n_sigma, Ptilde_n_sigma_fixed, Ptilde_n_sigma_estimate,
      Ptilde_sigma_fixed, Ptilde_sigma_estimate,
      Ptilde_scalarprior_mean,
      Ptilde_scalarprior_sd,
      Ptilde_raw_n_terms, Ptilde_raw_n_terms_fixed, Ptilde_raw_n_terms_estimate,
      Ptilde_raw_fixed, Ptilde_raw_estimate,
      Ptilde_re_start, Ptilde_re_end);
//   vector[n_geounit] tr_Ptilde = Ptilde_low + (1 - Ptilde_low)*inv_logit(get_mu(
// use probit
   vector[n_geounit] tr_Ptilde = //rep_vector(0.9999, n_geounit);
   Ptilde_low + (1 - Ptilde_low)*Phi_approx(get_mu(
      Ptilde_star,
      Ptilde_raw_n_terms,
      n_geounit,
      Ptilde_model_matrix_w, Ptilde_model_matrix_v, Ptilde_model_matrix_u));
  // omega
  vector[Omega_raw_n_terms] Omega_star = get_mu_star(
    Omega_n_sigma, Omega_n_sigma_fixed, Omega_n_sigma_estimate,
    Omega_sigma_fixed, Omega_sigma_estimate,
    Omega_scalarprior_mean, Omega_scalarprior_sd,
    Omega_raw_n_terms, Omega_raw_n_terms_fixed, Omega_raw_n_terms_estimate,
    Omega_raw_fixed, Omega_raw_estimate,
    Omega_re_start, Omega_re_end);
  vector[n_geounit] Omega = get_mu(
    Omega_star,
    Omega_raw_n_terms,
    n_geounit,
    Omega_model_matrix_w, Omega_model_matrix_v, Omega_model_matrix_u);

  // splines coeff
  matrix<lower=0>[Betas_n_sigma_estimate, k] Betas_sigma_estimate;
  if (Betas_n_sigma_estimate > 0) {
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 1] = reverse(Betas_sigma_estimate_reverse_1);
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 2] = reverse(Betas_sigma_estimate_reverse_2);
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 3] = reverse(Betas_sigma_estimate_reverse_3);
    Betas_sigma_estimate[1:Betas_n_sigma_estimate, 4] = reverse(Betas_sigma_estimate_reverse_4);
    //Betas_sigma_estimate[1:Betas_n_sigma_estimate, 5] = reverse(Betas_sigma_estimate_reverse_5);
    //   Betas_sigma_estimate[1:Betas_n_sigma_estimate, 6] = reverse(Betas_sigma_estimate_reverse_6);
    //   Betas_sigma_estimate[1:Betas_n_sigma_estimate, 7] = reverse(Betas_sigma_estimate_reverse_7);
  }
  matrix[Betas_raw_n_terms,Betas_k_terms] Betas_star = get_mudimhk_star(k ,
       Betas_n_sigma, Betas_n_sigma_fixed, Betas_n_sigma_estimate,
       Betas_sigma_fixed, Betas_sigma_estimate,
       Betas_scalarprior_mean, Betas_scalarprior_sd,
       Betas_raw_n_terms, Betas_raw_n_terms_fixed, Betas_raw_n_terms_estimate,
       Betas_raw_fixed, Betas_raw_estimate,
       Betas_re_start, Betas_re_end);
   matrix[n_geounit,Betas_k_terms] Betas = get_mudimhk(k,
       Betas_star,
       Betas_raw_n_terms,
       n_geounit,
       Betas_model_matrix_w, Betas_model_matrix_v, Betas_model_matrix_u);
  matrix[n_geounit,Betas_k_terms] tr_Betas_nonzero = Betas_lower_bound + (Betas_upper_bound - Betas_lower_bound) * inv_logit(Betas);


 // smoothing
  matrix[n_geounit, T] Epsilon;
  Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    Rho_fixed, Rho_estimate,
    Tau_fixed, Tau_estimate,
    Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_Eta_obs = process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing, Epsilon,
                        Omega, tr_Ptilde, tr_Betas_nonzero,
                        k, num_basis,   ext_knots,  spline_degree, add_shock);

    // asymptote
  vector[d_Ptilde_n_sigma_estimate] d_Ptilde_sigma_estimate = reverse(d_Ptilde_sigma_estimate_reverse);
  vector[d_Ptilde_raw_n_terms] d_Ptilde_star = get_mu_star(
      d_Ptilde_n_sigma, d_Ptilde_n_sigma_fixed, d_Ptilde_n_sigma_estimate,
      d_Ptilde_sigma_fixed, d_Ptilde_sigma_estimate,
      d_Ptilde_scalarprior_mean,
      d_Ptilde_scalarprior_sd,
      d_Ptilde_raw_n_terms, d_Ptilde_raw_n_terms_fixed, d_Ptilde_raw_n_terms_estimate,
      d_Ptilde_raw_fixed, d_Ptilde_raw_estimate,
      d_Ptilde_re_start, d_Ptilde_re_end);
//   vector[n_geounit] tr_d_Ptilde = d_Ptilde_low + (1 - d_Ptilde_low)*inv_logit(get_mu(
// use probit
   vector[n_geounit] tr_d_Ptilde = //rep_vector(0.9999, n_geounit);
   d_Ptilde_low + (1 - d_Ptilde_low)*Phi_approx(get_mu(
      d_Ptilde_star,
      d_Ptilde_raw_n_terms,
      n_geounit,
      d_Ptilde_model_matrix_w, d_Ptilde_model_matrix_v, d_Ptilde_model_matrix_u));
  // omega
  vector[d_Omega_raw_n_terms] d_Omega_star = get_mu_star(
    d_Omega_n_sigma, d_Omega_n_sigma_fixed, d_Omega_n_sigma_estimate,
    d_Omega_sigma_fixed, d_Omega_sigma_estimate,
    d_Omega_scalarprior_mean, d_Omega_scalarprior_sd,
    d_Omega_raw_n_terms, d_Omega_raw_n_terms_fixed, d_Omega_raw_n_terms_estimate,
    d_Omega_raw_fixed, d_Omega_raw_estimate,
    d_Omega_re_start, d_Omega_re_end);
  vector[n_geounit] d_Omega = get_mu(
    d_Omega_star,
    d_Omega_raw_n_terms,
    n_geounit,
    d_Omega_model_matrix_w, d_Omega_model_matrix_v, d_Omega_model_matrix_u);

  // splines coeff
  matrix<lower=0>[d_Betas_n_sigma_estimate, k] d_Betas_sigma_estimate;
  if (d_Betas_n_sigma_estimate > 0) {
    d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 1] = reverse(d_Betas_sigma_estimate_reverse_1);
    d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 2] = reverse(d_Betas_sigma_estimate_reverse_2);
    d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 3] = reverse(d_Betas_sigma_estimate_reverse_3);
    d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 4] = reverse(d_Betas_sigma_estimate_reverse_4);
    //d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 5] = reverse(d_Betas_sigma_estimate_reverse_5);
    //   d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 6] = reverse(d_Betas_sigma_estimate_reverse_6);
    //   d_Betas_sigma_estimate[1:d_Betas_n_sigma_estimate, 7] = reverse(d_Betas_sigma_estimate_reverse_7);
  }
  matrix[d_Betas_raw_n_terms,d_Betas_k_terms] d_Betas_star = get_mudimhk_star(k ,
       d_Betas_n_sigma, d_Betas_n_sigma_fixed, d_Betas_n_sigma_estimate,
       d_Betas_sigma_fixed, d_Betas_sigma_estimate,
       d_Betas_scalarprior_mean, d_Betas_scalarprior_sd,
       d_Betas_raw_n_terms, d_Betas_raw_n_terms_fixed, d_Betas_raw_n_terms_estimate,
       d_Betas_raw_fixed, d_Betas_raw_estimate,
       d_Betas_re_start, d_Betas_re_end);
   matrix[n_geounit,d_Betas_k_terms] d_Betas = get_mudimhk(k,
       d_Betas_star,
       d_Betas_raw_n_terms,
       n_geounit,
       d_Betas_model_matrix_w, d_Betas_model_matrix_v, d_Betas_model_matrix_u);
  matrix[n_geounit,d_Betas_k_terms] tr_d_Betas_nonzero = d_Betas_lower_bound + (d_Betas_upper_bound - d_Betas_lower_bound) * inv_logit(d_Betas);


 // smoothing
  matrix[n_geounit, T] d_Epsilon;
  d_Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    d_Rho_fixed, d_Rho_estimate,
    d_Tau_fixed, d_Tau_estimate,
    d_Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_d_Eta_obs = process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing, d_Epsilon,
                        d_Omega, tr_d_Ptilde, tr_d_Betas_nonzero,
                        k, num_basis,   ext_knots,  spline_degree, add_shock);

   

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


    // minor note that this includes obs that are NA, could subset instead
   array[add_dataoutliers * (1-DM2_fix_nonse)] real<lower=0> DM2_caux_dm_estimate = DM2_sqrt_caux_dm_estimate^2;
   vector[N] DM2_scale;
   if (add_dataoutliers){
    DM2_scale = get_scale(
        DM2_local_shrinkage_dm, DM2_nonse_fixed, DM2_global_shrinkage_dm_fixed, DM2_caux_dm_fixed,
        DM2_nonse_estimate, DM2_global_shrinkage_dm_estimate, DM2_caux_dm_estimate,
        N, DM2_s, isdhs, source, DM2_nooutlier, S, DM2_fix_nonse,
        DM2_any_bias, DM2_sdbias_fixed, DM2_sdbias_estimate);
   } else {
    DM2_scale = get_scale_nooutliers(
        DM2_nonse_fixed,
        DM2_nonse_estimate,
        N, DM2_s, isdhs, source, DM2_nooutlier, S, DM2_fix_nonse);
   }
     // pma
  real<lower = 0, upper = 1> DM2_rho_pma;
  if (DM2_fix_nonse) {
    DM2_rho_pma = DM2_rho_pma_fixed[1];
  } else {
    DM2_rho_pma = DM2_rho_pma_estimate[1];
  }
  // real<lower = 0, upper = 1> z_rho_pma;
  // if (add_z && fix_z) {
  //   z_rho_pma = z_rho_pma_fixed[1];
  // } else {
  //   z_rho_pma = z_rho_pma_estimate[1];
  // }


   

}

/////////////////////////////////////////////////////
model {

  // data model parameters, one set per dm
    DM1_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  DM1_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    DM1_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    DM1_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    DM1_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }

    DM2_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  DM2_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    DM2_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    DM2_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    DM2_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }

   

    // hierarchical parameters
  Ptilde_raw_estimate ~ std_normal();
  Omega_raw_estimate ~ std_normal();
  to_vector(Betas_raw_estimate) ~ std_normal();

  // variances
  //Ptilde_sigma_estimate ~ normal(0, Ptilde_prior_sd_sigma_estimate)T[0, positive_infinity()];
  Ptilde_sigma_estimate_reverse ~ normal(0, Ptilde_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  Omega_sigma_estimate ~ normal(0, Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(Betas_sigma_estimate) ~ normal(0, Betas_prior_sd_sigma_estimate)T[0, positive_infinity()];
   // one vector per spline coefficient
  to_vector(Betas_sigma_estimate_reverse_1) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(Betas_sigma_estimate_reverse_2) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(Betas_sigma_estimate_reverse_3) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(Betas_sigma_estimate_reverse_4) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(Betas_sigma_estimate_reverse_5) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(Betas_sigma_estimate_reverse_6) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(Betas_sigma_estimate_reverse_7) ~ normal(0, Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(Epsilon_innovation) ~ std_normal();
    if(fix_smoothing == 0) {
      Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }



    // hierarchical parameters
  d_Ptilde_raw_estimate ~ std_normal();
  d_Omega_raw_estimate ~ std_normal();
  to_vector(d_Betas_raw_estimate) ~ std_normal();

  // variances
  //d_Ptilde_sigma_estimate ~ normal(0, d_Ptilde_prior_sd_sigma_estimate)T[0, positive_infinity()];
  d_Ptilde_sigma_estimate_reverse ~ normal(0, d_Ptilde_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  d_Omega_sigma_estimate ~ normal(0, d_Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(d_Betas_sigma_estimate) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[0, positive_infinity()];
   // one vector per spline coefficient
  to_vector(d_Betas_sigma_estimate_reverse_1) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(d_Betas_sigma_estimate_reverse_2) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(d_Betas_sigma_estimate_reverse_3) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(d_Betas_sigma_estimate_reverse_4) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(d_Betas_sigma_estimate_reverse_5) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(d_Betas_sigma_estimate_reverse_6) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(d_Betas_sigma_estimate_reverse_7) ~ normal(0, d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(d_Epsilon_innovation) ~ std_normal();
    if(fix_smoothing == 0) {
      d_Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      d_Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }



   

  // using one dm block here for now to hardcode what mean is
    // model
  // fit to the data

  // get means for all, to use for pma and non-pma
  vector[N] modern;
  vector[N] logit_unmetovernonmodern;
  for (i in 1:N){
    real demand = inv_tr_eta(tr_d_Eta_obs[geo_unit[i], time[i]]);
    modern[i] = inv_tr_eta(tr_Eta_obs[geo_unit[i], time[i]]) * demand;
    logit_unmetovernonmodern[i] = logit((demand - modern[i])/(1-modern[i]));
  }

  // non-pma data
  for(i in 1:N) {
    if(held_out[i] == 0 && ispma[i] == 0) {
      if(DM1_obs_isna[i] == 0) {
          DM1_y[i] ~ normal(logit(modern[i]), DM1_scale[i]);
      }
      if(DM2_obs_isna[i] == 0) {
         DM2_y[i] ~ normal(logit_unmetovernonmodern[i], DM2_scale[i]);
      }
    }
  }
  // pma
  // no densities, constrained between 0 and 1
  // not yet for held out (held out just filtered out)
  // not yet for agregates
  // same for unmet
  for (index_country in pma_country_indices){
    array[npma_c[index_country]] int indices_inc = indices_pma[nstart_pma_c[index_country]:
          (nstart_pma_c[index_country]+npma_c[index_country] - 1)];
    matrix[npma_c[index_country], npma_c[index_country]] DM1_cov_matrix_pma;
    matrix[npma_c[index_country], npma_c[index_country]] DM2_cov_matrix_pma;
    for (j1 in 1:npma_c[index_country]){
      for (j2 in 1:npma_c[index_country]){
        int j1_index = indices_inc[j1];
        int j2_index = indices_inc[j2];
        DM1_cov_matrix_pma[j1, j2] =
            DM1_scale[j1_index]* DM1_scale[j2_index] *
              DM1_rho_pma^(abs(exact_time[j1_index] - exact_time[j2_index]));
        DM2_cov_matrix_pma[j1, j2] =
            DM2_scale[j1_index]* DM2_scale[j2_index] *
              DM2_rho_pma^(abs(exact_time[j1_index] - exact_time[j2_index]));
    }}
    DM1_y[indices_inc] ~ multi_normal(logit(modern[indices_inc]), DM1_cov_matrix_pma);
    DM2_y[indices_inc] ~ multi_normal(logit_unmetovernonmodern[indices_inc], DM2_cov_matrix_pma);
  }






  // subnational
   if (!fix_subnat_corr){
    rho_correlationeps_estimate[1] ~ uniform(0, 1);
    }

     

}


generated quantities {

    matrix[n_geounit, T] Eta = inv_tr_eta_matrix(
            process_model_outsideobs(tr_Eta_obs,
                      n_geounit, T,  t_min, t_max,
                      smoothing, Epsilon,
                      tr_Ptilde, tr_Betas_nonzero,
                      k, num_basis,   ext_knots,  spline_degree, add_shock));

    matrix[n_geounit, T] d_Eta = inv_tr_eta_matrix(
            process_model_outsideobs(tr_d_Eta_obs,
                      n_geounit, T,  t_min, t_max,
                      smoothing, d_Epsilon,
                      tr_d_Ptilde, tr_d_Betas_nonzero,
                      k, num_basis,   ext_knots,  spline_degree, add_shock));

   
   
   
   

}






