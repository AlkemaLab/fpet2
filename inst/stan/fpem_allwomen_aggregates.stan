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

// for unmarried, add unmarried_
// note that geounit_index need to match between married and unmarried
  int unmarried_N; // Number of observations
  int unmarried_S; // Number of sources
  array[unmarried_N] int<lower=1, upper=T> unmarried_time;     // Time of each observation
  array[unmarried_N] int<lower=0, upper=n_geounit> unmarried_geo_unit; // Geographic unit of each observation, 0 if observation is for an aggregated unit
  // updated to lower of 0 for subnational runs
  array[unmarried_N] int<lower=0, upper=S> unmarried_source;   // Source of each observation
  array[unmarried_N] int<lower=0, upper=1> unmarried_held_out; // Whether to hold out each observation
  array[unmarried_N] int<lower=0> unmarried_isdhs;            // DHS?
  // pma data
  array[unmarried_N] real<lower=0> unmarried_ispma;            // PMA?
  array[unmarried_N] real<lower=0> unmarried_exact_time;
  int<lower=0> unmarried_N_pma_country_indices;
  array[unmarried_N_pma_country_indices] int<lower=1, upper = n_geounit> unmarried_pma_country_indices;
  array[n_geounit] int unmarried_npma_c;
  array[n_geounit] int unmarried_nstart_pma_c;
  int<lower=0> unmarried_N_pma_obs_indices;
  array[unmarried_N_pma_obs_indices] int unmarried_indices_pma;

  matrix[n_geounit, T] prop_married_rt; // proportion married for each region-year
  vector[T] prop_married_t; // at national level, same as prop_married_rt for national run



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


  real unmarried_Betas_lower_bound;
  real unmarried_Betas_upper_bound;
  real unmarried_Ptilde_low;

  int<lower=0> unmarried_Ptilde_raw_n_terms;
  int<lower=0> unmarried_Ptilde_raw_n_terms_fixed;
  int<lower=0> unmarried_Ptilde_raw_n_terms_estimate;
  int<lower=0> unmarried_Ptilde_n_sigma;
  int<lower=0> unmarried_Ptilde_n_sigma_fixed;
  int<lower=0> unmarried_Ptilde_n_sigma_estimate;
  array[unmarried_Ptilde_n_sigma + 1] int<lower=1, upper=unmarried_Ptilde_raw_n_terms> unmarried_Ptilde_re_start;
  array[unmarried_Ptilde_n_sigma + 1] int<lower=1, upper=unmarried_Ptilde_raw_n_terms> unmarried_Ptilde_re_end;
  matrix[n_geounit, unmarried_Ptilde_raw_n_terms] unmarried_Ptilde_model_matrix;
  real<lower = verysmallnumber> unmarried_Ptilde_scalarprior_sd;
  real unmarried_Ptilde_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_Ptilde_prior_sd_sigma_estimate;
  vector[unmarried_Ptilde_raw_n_terms_fixed] unmarried_Ptilde_raw_fixed;
  vector<lower=0>[unmarried_Ptilde_n_sigma_fixed] unmarried_Ptilde_sigma_fixed;

  int<lower=0> unmarried_Omega_raw_n_terms;
  int<lower=0> unmarried_Omega_n_sigma;
  int<lower=0> unmarried_Omega_raw_n_terms_fixed;
  int<lower=0> unmarried_Omega_raw_n_terms_estimate;
  int<lower=0> unmarried_Omega_n_sigma_fixed;
  int<lower=0> unmarried_Omega_n_sigma_estimate;
  array[unmarried_Omega_n_sigma + 1] int<lower=1, upper=unmarried_Omega_raw_n_terms> unmarried_Omega_re_start;
  array[unmarried_Omega_n_sigma + 1] int<lower=1, upper=unmarried_Omega_raw_n_terms> unmarried_Omega_re_end;
  matrix[n_geounit, unmarried_Omega_raw_n_terms] unmarried_Omega_model_matrix;
  real<lower = verysmallnumber> unmarried_Omega_scalarprior_sd;
  real unmarried_Omega_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_Omega_prior_sd_sigma_estimate;
  vector[unmarried_Omega_raw_n_terms_fixed] unmarried_Omega_raw_fixed;
  vector<lower=0>[unmarried_Omega_n_sigma_fixed] unmarried_Omega_sigma_fixed;

  int unmarried_Betas_k_terms; // same as k; to do clean to have just 1
  int<lower=0> unmarried_Betas_raw_n_terms;
  int<lower=0> unmarried_Betas_raw_n_terms_fixed;
  int<lower=0> unmarried_Betas_raw_n_terms_estimate;
  int<lower=0> unmarried_Betas_n_sigma;
  int<lower=0> unmarried_Betas_n_sigma_fixed;
  int<lower=0> unmarried_Betas_n_sigma_estimate;
  array[unmarried_Betas_n_sigma + 1] int<lower=1, upper = unmarried_Betas_raw_n_terms> unmarried_Betas_re_start;
  array[unmarried_Betas_n_sigma + 1] int<lower=1, upper = unmarried_Betas_raw_n_terms> unmarried_Betas_re_end;
  real<lower = verysmallnumber> unmarried_Betas_scalarprior_sd;
  real unmarried_Betas_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_Betas_prior_sd_sigma_estimate;
  matrix[n_geounit, unmarried_Betas_raw_n_terms] unmarried_Betas_model_matrix;
  matrix[unmarried_Betas_raw_n_terms_fixed, k] unmarried_Betas_raw_fixed;
  matrix<lower=0>[unmarried_Betas_n_sigma_fixed, k] unmarried_Betas_sigma_fixed;


    // // these can be used as max's when adding a level
  // // not yet used
  // real unmarried_Ptilde_sigma_max;
  // real unmarried_Omega_sigma_max;
  // real a_sigma_max_1;
  // real a_sigma_max_2;
  // real a_sigma_max_3;
  // real a_sigma_max_4;


    // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[fix_smoothing] real<lower=0, upper=1> unmarried_Rho_fixed;
  array[fix_smoothing] real<lower=0> unmarried_Tau_fixed;


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


  real unmarried_d_Betas_lower_bound;
  real unmarried_d_Betas_upper_bound;
  real unmarried_d_Ptilde_low;

  int<lower=0> unmarried_d_Ptilde_raw_n_terms;
  int<lower=0> unmarried_d_Ptilde_raw_n_terms_fixed;
  int<lower=0> unmarried_d_Ptilde_raw_n_terms_estimate;
  int<lower=0> unmarried_d_Ptilde_n_sigma;
  int<lower=0> unmarried_d_Ptilde_n_sigma_fixed;
  int<lower=0> unmarried_d_Ptilde_n_sigma_estimate;
  array[unmarried_d_Ptilde_n_sigma + 1] int<lower=1, upper=unmarried_d_Ptilde_raw_n_terms> unmarried_d_Ptilde_re_start;
  array[unmarried_d_Ptilde_n_sigma + 1] int<lower=1, upper=unmarried_d_Ptilde_raw_n_terms> unmarried_d_Ptilde_re_end;
  matrix[n_geounit, unmarried_d_Ptilde_raw_n_terms] unmarried_d_Ptilde_model_matrix;
  real<lower = verysmallnumber> unmarried_d_Ptilde_scalarprior_sd;
  real unmarried_d_Ptilde_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_d_Ptilde_prior_sd_sigma_estimate;
  vector[unmarried_d_Ptilde_raw_n_terms_fixed] unmarried_d_Ptilde_raw_fixed;
  vector<lower=0>[unmarried_d_Ptilde_n_sigma_fixed] unmarried_d_Ptilde_sigma_fixed;

  int<lower=0> unmarried_d_Omega_raw_n_terms;
  int<lower=0> unmarried_d_Omega_n_sigma;
  int<lower=0> unmarried_d_Omega_raw_n_terms_fixed;
  int<lower=0> unmarried_d_Omega_raw_n_terms_estimate;
  int<lower=0> unmarried_d_Omega_n_sigma_fixed;
  int<lower=0> unmarried_d_Omega_n_sigma_estimate;
  array[unmarried_d_Omega_n_sigma + 1] int<lower=1, upper=unmarried_d_Omega_raw_n_terms> unmarried_d_Omega_re_start;
  array[unmarried_d_Omega_n_sigma + 1] int<lower=1, upper=unmarried_d_Omega_raw_n_terms> unmarried_d_Omega_re_end;
  matrix[n_geounit, unmarried_d_Omega_raw_n_terms] unmarried_d_Omega_model_matrix;
  real<lower = verysmallnumber> unmarried_d_Omega_scalarprior_sd;
  real unmarried_d_Omega_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_d_Omega_prior_sd_sigma_estimate;
  vector[unmarried_d_Omega_raw_n_terms_fixed] unmarried_d_Omega_raw_fixed;
  vector<lower=0>[unmarried_d_Omega_n_sigma_fixed] unmarried_d_Omega_sigma_fixed;

  int unmarried_d_Betas_k_terms; // same as k; to do clean to have just 1
  int<lower=0> unmarried_d_Betas_raw_n_terms;
  int<lower=0> unmarried_d_Betas_raw_n_terms_fixed;
  int<lower=0> unmarried_d_Betas_raw_n_terms_estimate;
  int<lower=0> unmarried_d_Betas_n_sigma;
  int<lower=0> unmarried_d_Betas_n_sigma_fixed;
  int<lower=0> unmarried_d_Betas_n_sigma_estimate;
  array[unmarried_d_Betas_n_sigma + 1] int<lower=1, upper = unmarried_d_Betas_raw_n_terms> unmarried_d_Betas_re_start;
  array[unmarried_d_Betas_n_sigma + 1] int<lower=1, upper = unmarried_d_Betas_raw_n_terms> unmarried_d_Betas_re_end;
  real<lower = verysmallnumber> unmarried_d_Betas_scalarprior_sd;
  real unmarried_d_Betas_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_d_Betas_prior_sd_sigma_estimate;
  matrix[n_geounit, unmarried_d_Betas_raw_n_terms] unmarried_d_Betas_model_matrix;
  matrix[unmarried_d_Betas_raw_n_terms_fixed, k] unmarried_d_Betas_raw_fixed;
  matrix<lower=0>[unmarried_d_Betas_n_sigma_fixed, k] unmarried_d_Betas_sigma_fixed;


    // // these can be used as max's when adding a level
  // // not yet used
  // real unmarried_d_Ptilde_sigma_max;
  // real unmarried_d_Omega_sigma_max;
  // real a_sigma_max_1;
  // real a_sigma_max_2;
  // real a_sigma_max_3;
  // real a_sigma_max_4;


    // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[fix_smoothing] real<lower=0, upper=1> unmarried_d_Rho_fixed;
  array[fix_smoothing] real<lower=0> unmarried_d_Tau_fixed;


  
  int<lower=0, upper=1> z_fix_smoothing;

  int<lower=0> z_Omega_raw_n_terms;
  int<lower=0> z_Omega_n_sigma;
  int<lower=0> z_Omega_raw_n_terms_fixed;
  int<lower=0> z_Omega_raw_n_terms_estimate;
  int<lower=0> z_Omega_n_sigma_fixed;
  int<lower=0> z_Omega_n_sigma_estimate;
  array[z_Omega_n_sigma + 1] int<lower=1, upper=z_Omega_raw_n_terms> z_Omega_re_start;
  array[z_Omega_n_sigma + 1] int<lower=1, upper=z_Omega_raw_n_terms> z_Omega_re_end;
  matrix[n_geounit, z_Omega_raw_n_terms] z_Omega_model_matrix;
  real<lower = verysmallnumber> z_Omega_scalarprior_sd;
  real z_Omega_scalarprior_mean;
  real<lower = verysmallnumber> z_Omega_prior_sd_sigma_estimate;
  vector[z_Omega_raw_n_terms_fixed] z_Omega_raw_fixed;
  vector<lower=0>[z_Omega_n_sigma_fixed] z_Omega_sigma_fixed;


  // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[z_fix_smoothing] real<lower=0, upper=1> z_Rho_fixed;
  array[z_fix_smoothing] real<lower=0> z_Tau_fixed;



  int<lower=0, upper=1> unmarried_z_fix_smoothing;

  int<lower=0> unmarried_z_Omega_raw_n_terms;
  int<lower=0> unmarried_z_Omega_n_sigma;
  int<lower=0> unmarried_z_Omega_raw_n_terms_fixed;
  int<lower=0> unmarried_z_Omega_raw_n_terms_estimate;
  int<lower=0> unmarried_z_Omega_n_sigma_fixed;
  int<lower=0> unmarried_z_Omega_n_sigma_estimate;
  array[unmarried_z_Omega_n_sigma + 1] int<lower=1, upper=unmarried_z_Omega_raw_n_terms> unmarried_z_Omega_re_start;
  array[unmarried_z_Omega_n_sigma + 1] int<lower=1, upper=unmarried_z_Omega_raw_n_terms> unmarried_z_Omega_re_end;
  matrix[n_geounit, unmarried_z_Omega_raw_n_terms] unmarried_z_Omega_model_matrix;
  real<lower = verysmallnumber> unmarried_z_Omega_scalarprior_sd;
  real unmarried_z_Omega_scalarprior_mean;
  real<lower = verysmallnumber> unmarried_z_Omega_prior_sd_sigma_estimate;
  vector[unmarried_z_Omega_raw_n_terms_fixed] unmarried_z_Omega_raw_fixed;
  vector<lower=0>[unmarried_z_Omega_n_sigma_fixed] unmarried_z_Omega_sigma_fixed;


  // only relevant if smoothing = 1, should be 0 if smoothing is 0
  array[unmarried_z_fix_smoothing] real<lower=0, upper=1> unmarried_z_Rho_fixed;
  array[unmarried_z_fix_smoothing] real<lower=0> unmarried_z_Tau_fixed;





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
  int<lower=0, upper=1> unmarried_DM1_fix_nonse; // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  vector[unmarried_N] unmarried_DM1_y;                // Observations
  vector[unmarried_N] unmarried_DM1_s;                // unmarried_Standard deviation
  array[unmarried_N] int<lower=0> unmarried_DM1_obs_isna;
  array[unmarried_N] int<lower=0> unmarried_DM1_any_bias;
  array[unmarried_N] int<lower=0> unmarried_DM1_nooutlier;                // should outlier be included?
  array[unmarried_DM1_fix_nonse ? unmarried_S : 0] real<lower=0> unmarried_DM1_nonse_fixed;
  array[add_dataoutliers*unmarried_DM1_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM1_global_shrinkage_dm_fixed;
  array[add_dataoutliers*unmarried_DM1_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM1_caux_dm_fixed;
  array[unmarried_DM1_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM1_sdbias_fixed;
  array[unmarried_DM1_fix_nonse ? 1 : 0] real<lower=0, upper = 1> unmarried_DM1_rho_pma_fixed;
  //array[unmarried_DM1_fix_nonse*add_z*fix_z ? 1 : 0] real<lower=0, upper = 1> z_rho_pma_fixed;

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

  // one set per dm
  int<lower=0, upper=1> unmarried_DM2_fix_nonse; // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  vector[unmarried_N] unmarried_DM2_y;                // Observations
  vector[unmarried_N] unmarried_DM2_s;                // unmarried_Standard deviation
  array[unmarried_N] int<lower=0> unmarried_DM2_obs_isna;
  array[unmarried_N] int<lower=0> unmarried_DM2_any_bias;
  array[unmarried_N] int<lower=0> unmarried_DM2_nooutlier;                // should outlier be included?
  array[unmarried_DM2_fix_nonse ? unmarried_S : 0] real<lower=0> unmarried_DM2_nonse_fixed;
  array[add_dataoutliers*unmarried_DM2_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM2_global_shrinkage_dm_fixed;
  array[add_dataoutliers*unmarried_DM2_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM2_caux_dm_fixed;
  array[unmarried_DM2_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM2_sdbias_fixed;
  array[unmarried_DM2_fix_nonse ? 1 : 0] real<lower=0, upper = 1> unmarried_DM2_rho_pma_fixed;
  //array[unmarried_DM2_fix_nonse*add_z*fix_z ? 1 : 0] real<lower=0, upper = 1> z_rho_pma_fixed;

    // one set per dm
  int<lower=0, upper=1> DM3_fix_nonse; // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  vector[N] DM3_y;                // Observations
  vector[N] DM3_s;                // Standard deviation
  array[N] int<lower=0> DM3_obs_isna;
  array[N] int<lower=0> DM3_any_bias;
  array[N] int<lower=0> DM3_nooutlier;                // should outlier be included?
  array[DM3_fix_nonse ? S : 0] real<lower=0> DM3_nonse_fixed;
  array[add_dataoutliers*DM3_fix_nonse ? 1 : 0] real<lower=0> DM3_global_shrinkage_dm_fixed;
  array[add_dataoutliers*DM3_fix_nonse ? 1 : 0] real<lower=0> DM3_caux_dm_fixed;
  array[DM3_fix_nonse ? 1 : 0] real<lower=0> DM3_sdbias_fixed;
  array[DM3_fix_nonse ? 1 : 0] real<lower=0, upper = 1> DM3_rho_pma_fixed;
  //array[DM3_fix_nonse*add_z*fix_z ? 1 : 0] real<lower=0, upper = 1> z_rho_pma_fixed;

  // one set per dm
  int<lower=0, upper=1> unmarried_DM3_fix_nonse; // 0: estimate nonse. 1: fix nonse at estimates from previous fit.
  vector[unmarried_N] unmarried_DM3_y;                // Observations
  vector[unmarried_N] unmarried_DM3_s;                // unmarried_Standard deviation
  array[unmarried_N] int<lower=0> unmarried_DM3_obs_isna;
  array[unmarried_N] int<lower=0> unmarried_DM3_any_bias;
  array[unmarried_N] int<lower=0> unmarried_DM3_nooutlier;                // should outlier be included?
  array[unmarried_DM3_fix_nonse ? unmarried_S : 0] real<lower=0> unmarried_DM3_nonse_fixed;
  array[add_dataoutliers*unmarried_DM3_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM3_global_shrinkage_dm_fixed;
  array[add_dataoutliers*unmarried_DM3_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM3_caux_dm_fixed;
  array[unmarried_DM3_fix_nonse ? 1 : 0] real<lower=0> unmarried_DM3_sdbias_fixed;
  array[unmarried_DM3_fix_nonse ? 1 : 0] real<lower=0, upper = 1> unmarried_DM3_rho_pma_fixed;
  //array[unmarried_DM3_fix_nonse*add_z*fix_z ? 1 : 0] real<lower=0, upper = 1> z_rho_pma_fixed;


    // subnational aggregates
  matrix[T, n_geounit] geo_unit_natpop_weight_tr; // Weights used to aggregate geo units for national aggregates. Each column should sum to 1

  // subnational aggregates
  matrix[T, n_geounit] unmarried_geo_unit_natpop_weight_tr; // Weights used to aggregate geo units for national aggregates. Each column should sum to 1


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

 vector[rows(csr_extract_w(unmarried_Ptilde_model_matrix))] unmarried_Ptilde_model_matrix_w    = csr_extract_w(unmarried_Ptilde_model_matrix);
  array[size(csr_extract_v(unmarried_Ptilde_model_matrix))] int unmarried_Ptilde_model_matrix_v = csr_extract_v(unmarried_Ptilde_model_matrix);
  array[size(csr_extract_u(unmarried_Ptilde_model_matrix))] int unmarried_Ptilde_model_matrix_u = csr_extract_u(unmarried_Ptilde_model_matrix);

  vector[rows(csr_extract_w(unmarried_Omega_model_matrix))] unmarried_Omega_model_matrix_w     = csr_extract_w(unmarried_Omega_model_matrix);
  array[size(csr_extract_v(unmarried_Omega_model_matrix))] int unmarried_Omega_model_matrix_v  = csr_extract_v(unmarried_Omega_model_matrix);
  array[size(csr_extract_u(unmarried_Omega_model_matrix))] int unmarried_Omega_model_matrix_u  = csr_extract_u(unmarried_Omega_model_matrix);

  vector[rows(csr_extract_w(unmarried_Betas_model_matrix))] unmarried_Betas_model_matrix_w    = csr_extract_w(unmarried_Betas_model_matrix);
  array[size(csr_extract_v(unmarried_Betas_model_matrix))] int unmarried_Betas_model_matrix_v = csr_extract_v(unmarried_Betas_model_matrix);
  array[size(csr_extract_u(unmarried_Betas_model_matrix))] int unmarried_Betas_model_matrix_u = csr_extract_u(unmarried_Betas_model_matrix);

   vector[rows(csr_extract_w(d_Ptilde_model_matrix))] d_Ptilde_model_matrix_w    = csr_extract_w(d_Ptilde_model_matrix);
  array[size(csr_extract_v(d_Ptilde_model_matrix))] int d_Ptilde_model_matrix_v = csr_extract_v(d_Ptilde_model_matrix);
  array[size(csr_extract_u(d_Ptilde_model_matrix))] int d_Ptilde_model_matrix_u = csr_extract_u(d_Ptilde_model_matrix);

  vector[rows(csr_extract_w(d_Omega_model_matrix))] d_Omega_model_matrix_w     = csr_extract_w(d_Omega_model_matrix);
  array[size(csr_extract_v(d_Omega_model_matrix))] int d_Omega_model_matrix_v  = csr_extract_v(d_Omega_model_matrix);
  array[size(csr_extract_u(d_Omega_model_matrix))] int d_Omega_model_matrix_u  = csr_extract_u(d_Omega_model_matrix);

  vector[rows(csr_extract_w(d_Betas_model_matrix))] d_Betas_model_matrix_w    = csr_extract_w(d_Betas_model_matrix);
  array[size(csr_extract_v(d_Betas_model_matrix))] int d_Betas_model_matrix_v = csr_extract_v(d_Betas_model_matrix);
  array[size(csr_extract_u(d_Betas_model_matrix))] int d_Betas_model_matrix_u = csr_extract_u(d_Betas_model_matrix);

 vector[rows(csr_extract_w(unmarried_d_Ptilde_model_matrix))] unmarried_d_Ptilde_model_matrix_w    = csr_extract_w(unmarried_d_Ptilde_model_matrix);
  array[size(csr_extract_v(unmarried_d_Ptilde_model_matrix))] int unmarried_d_Ptilde_model_matrix_v = csr_extract_v(unmarried_d_Ptilde_model_matrix);
  array[size(csr_extract_u(unmarried_d_Ptilde_model_matrix))] int unmarried_d_Ptilde_model_matrix_u = csr_extract_u(unmarried_d_Ptilde_model_matrix);

  vector[rows(csr_extract_w(unmarried_d_Omega_model_matrix))] unmarried_d_Omega_model_matrix_w     = csr_extract_w(unmarried_d_Omega_model_matrix);
  array[size(csr_extract_v(unmarried_d_Omega_model_matrix))] int unmarried_d_Omega_model_matrix_v  = csr_extract_v(unmarried_d_Omega_model_matrix);
  array[size(csr_extract_u(unmarried_d_Omega_model_matrix))] int unmarried_d_Omega_model_matrix_u  = csr_extract_u(unmarried_d_Omega_model_matrix);

  vector[rows(csr_extract_w(unmarried_d_Betas_model_matrix))] unmarried_d_Betas_model_matrix_w    = csr_extract_w(unmarried_d_Betas_model_matrix);
  array[size(csr_extract_v(unmarried_d_Betas_model_matrix))] int unmarried_d_Betas_model_matrix_v = csr_extract_v(unmarried_d_Betas_model_matrix);
  array[size(csr_extract_u(unmarried_d_Betas_model_matrix))] int unmarried_d_Betas_model_matrix_u = csr_extract_u(unmarried_d_Betas_model_matrix);

  
  vector[rows(csr_extract_w(z_Omega_model_matrix))] z_Omega_model_matrix_w     = csr_extract_w(z_Omega_model_matrix);
  array[size(csr_extract_v(z_Omega_model_matrix))] int z_Omega_model_matrix_v  = csr_extract_v(z_Omega_model_matrix);
  array[size(csr_extract_u(z_Omega_model_matrix))] int z_Omega_model_matrix_u  = csr_extract_u(z_Omega_model_matrix);



  vector[rows(csr_extract_w(unmarried_z_Omega_model_matrix))] unmarried_z_Omega_model_matrix_w     = csr_extract_w(unmarried_z_Omega_model_matrix);
  array[size(csr_extract_v(unmarried_z_Omega_model_matrix))] int unmarried_z_Omega_model_matrix_v  = csr_extract_v(unmarried_z_Omega_model_matrix);
  array[size(csr_extract_u(unmarried_z_Omega_model_matrix))] int unmarried_z_Omega_model_matrix_u  = csr_extract_u(unmarried_z_Omega_model_matrix);




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
  matrix[unmarried_Betas_raw_n_terms_estimate, k] unmarried_Betas_raw_estimate;
  // matrix<lower = 0>[a_n_sigma_estimate, k>3 ? 3 : k] a_sigma_estimate;
  // with ordered variances
  positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_1;
  positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_2;
  positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_3;
  positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_4;
  //positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_5;
  // positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_6;
  // positive_ordered [unmarried_Betas_n_sigma_estimate] unmarried_Betas_sigma_estimate_reverse_7;

  vector[unmarried_Ptilde_raw_n_terms_estimate] unmarried_Ptilde_raw_estimate;
  //vector<lower=verysmallnumber>[unmarried_Ptilde_n_sigma_estimate] unmarried_Ptilde_sigma_estimate;
  positive_ordered [unmarried_Ptilde_n_sigma_estimate] unmarried_Ptilde_sigma_estimate_reverse;

  vector[unmarried_Omega_raw_n_terms_estimate] unmarried_Omega_raw_estimate;
  vector<lower=0>[unmarried_Omega_n_sigma_estimate] unmarried_Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] unmarried_Epsilon_innovation;
  array[smoothing * (1 - fix_smoothing)] real<lower=0, upper=1> unmarried_Rho_estimate;
  array[smoothing * (1 - fix_smoothing)] real<lower=0> unmarried_Tau_estimate;

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

  // splines coefficients
  matrix[unmarried_d_Betas_raw_n_terms_estimate, k] unmarried_d_Betas_raw_estimate;
  // matrix<lower = 0>[a_n_sigma_estimate, k>3 ? 3 : k] a_sigma_estimate;
  // with ordered variances
  positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_1;
  positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_2;
  positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_3;
  positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_4;
  //positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_5;
  // positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_6;
  // positive_ordered [unmarried_d_Betas_n_sigma_estimate] unmarried_d_Betas_sigma_estimate_reverse_7;

  vector[unmarried_d_Ptilde_raw_n_terms_estimate] unmarried_d_Ptilde_raw_estimate;
  //vector<lower=verysmallnumber>[unmarried_d_Ptilde_n_sigma_estimate] unmarried_d_Ptilde_sigma_estimate;
  positive_ordered [unmarried_d_Ptilde_n_sigma_estimate] unmarried_d_Ptilde_sigma_estimate_reverse;

  vector[unmarried_d_Omega_raw_n_terms_estimate] unmarried_d_Omega_raw_estimate;
  vector<lower=0>[unmarried_d_Omega_n_sigma_estimate] unmarried_d_Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] unmarried_d_Epsilon_innovation;
  array[smoothing * (1 - fix_smoothing)] real<lower=0, upper=1> unmarried_d_Rho_estimate;
  array[smoothing * (1 - fix_smoothing)] real<lower=0> unmarried_d_Tau_estimate;

   
  vector[z_Omega_raw_n_terms_estimate] z_Omega_raw_estimate;
  vector<lower=0>[z_Omega_n_sigma_estimate] z_Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] z_Epsilon_innovation;
  array[smoothing * (1 - z_fix_smoothing)] real<lower=0, upper=1> z_Rho_estimate;
  array[smoothing * (1 - z_fix_smoothing)] real<lower=0> z_Tau_estimate;


  vector[unmarried_z_Omega_raw_n_terms_estimate] unmarried_z_Omega_raw_estimate;
  vector<lower=0>[unmarried_z_Omega_n_sigma_estimate] unmarried_z_Omega_sigma_estimate;

  matrix[n_geounit * smoothing, T * smoothing] unmarried_z_Epsilon_innovation;
  array[smoothing * (1 - unmarried_z_fix_smoothing)] real<lower=0, upper=1> unmarried_z_Rho_estimate;
  array[smoothing * (1 - unmarried_z_fix_smoothing)] real<lower=0> unmarried_z_Tau_estimate;


  // Data model, one set per dm
    array[DM1_fix_nonse ? 0 : S] real<lower=0> DM1_nonse_estimate;
  array[DM1_fix_nonse ? 0 : 1] real<lower=0> DM1_sdbias_estimate;
  array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-DM1_fix_nonse)] real<lower=0> DM1_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*N] DM1_local_shrinkage_dm;
  array[(1-DM1_fix_nonse)] real<lower = 0, upper = 0.99> DM1_rho_pma_estimate;
 // array[(1-DM1_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


  array[unmarried_DM1_fix_nonse ? 0 : unmarried_S] real<lower=0> unmarried_DM1_nonse_estimate;
  array[unmarried_DM1_fix_nonse ? 0 : 1] real<lower=0> unmarried_DM1_sdbias_estimate;
  array[add_dataoutliers * (1-unmarried_DM1_fix_nonse)] real<lower=0> unmarried_DM1_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-unmarried_DM1_fix_nonse)] real<lower=0> unmarried_DM1_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*unmarried_N] unmarried_DM1_local_shrinkage_dm;
  array[(1-unmarried_DM1_fix_nonse)] real<lower = 0, upper = 0.99> unmarried_DM1_rho_pma_estimate;
 // array[(1-unmarried_DM1_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


    array[DM2_fix_nonse ? 0 : S] real<lower=0> DM2_nonse_estimate;
  array[DM2_fix_nonse ? 0 : 1] real<lower=0> DM2_sdbias_estimate;
  array[add_dataoutliers * (1-DM2_fix_nonse)] real<lower=0> DM2_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-DM2_fix_nonse)] real<lower=0> DM2_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*N] DM2_local_shrinkage_dm;
  array[(1-DM2_fix_nonse)] real<lower = 0, upper = 0.99> DM2_rho_pma_estimate;
 // array[(1-DM2_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


  array[unmarried_DM2_fix_nonse ? 0 : unmarried_S] real<lower=0> unmarried_DM2_nonse_estimate;
  array[unmarried_DM2_fix_nonse ? 0 : 1] real<lower=0> unmarried_DM2_sdbias_estimate;
  array[add_dataoutliers * (1-unmarried_DM2_fix_nonse)] real<lower=0> unmarried_DM2_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-unmarried_DM2_fix_nonse)] real<lower=0> unmarried_DM2_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*unmarried_N] unmarried_DM2_local_shrinkage_dm;
  array[(1-unmarried_DM2_fix_nonse)] real<lower = 0, upper = 0.99> unmarried_DM2_rho_pma_estimate;
 // array[(1-unmarried_DM2_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


     array[DM3_fix_nonse ? 0 : S] real<lower=0> DM3_nonse_estimate;
  array[DM3_fix_nonse ? 0 : 1] real<lower=0> DM3_sdbias_estimate;
  array[add_dataoutliers * (1-DM3_fix_nonse)] real<lower=0> DM3_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-DM3_fix_nonse)] real<lower=0> DM3_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*N] DM3_local_shrinkage_dm;
  array[(1-DM3_fix_nonse)] real<lower = 0, upper = 0.99> DM3_rho_pma_estimate;
 // array[(1-DM3_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;


  array[unmarried_DM3_fix_nonse ? 0 : unmarried_S] real<lower=0> unmarried_DM3_nonse_estimate;
  array[unmarried_DM3_fix_nonse ? 0 : 1] real<lower=0> unmarried_DM3_sdbias_estimate;
  array[add_dataoutliers * (1-unmarried_DM3_fix_nonse)] real<lower=0> unmarried_DM3_global_shrinkage_dm_estimate;
  array[add_dataoutliers * (1-unmarried_DM3_fix_nonse)] real<lower=0> unmarried_DM3_sqrt_caux_dm_estimate;
  vector<lower=0>[add_dataoutliers*unmarried_N] unmarried_DM3_local_shrinkage_dm;
  array[(1-unmarried_DM3_fix_nonse)] real<lower = 0, upper = 0.99> unmarried_DM3_rho_pma_estimate;
 // array[(1-unmarried_DM3_fix_nonse*fix_z)] real<lower = 0, upper = 0.99> z_rho_pma_estimate;



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
  vector[unmarried_Ptilde_n_sigma_estimate] unmarried_Ptilde_sigma_estimate = reverse(unmarried_Ptilde_sigma_estimate_reverse);
  vector[unmarried_Ptilde_raw_n_terms] unmarried_Ptilde_star = get_mu_star(
      unmarried_Ptilde_n_sigma, unmarried_Ptilde_n_sigma_fixed, unmarried_Ptilde_n_sigma_estimate,
      unmarried_Ptilde_sigma_fixed, unmarried_Ptilde_sigma_estimate,
      unmarried_Ptilde_scalarprior_mean,
      unmarried_Ptilde_scalarprior_sd,
      unmarried_Ptilde_raw_n_terms, unmarried_Ptilde_raw_n_terms_fixed, unmarried_Ptilde_raw_n_terms_estimate,
      unmarried_Ptilde_raw_fixed, unmarried_Ptilde_raw_estimate,
      unmarried_Ptilde_re_start, unmarried_Ptilde_re_end);
//   vector[n_geounit] tr_unmarried_Ptilde = unmarried_Ptilde_low + (1 - unmarried_Ptilde_low)*inv_logit(get_mu(
// use probit
   vector[n_geounit] tr_unmarried_Ptilde = //rep_vector(0.9999, n_geounit);
   unmarried_Ptilde_low + (1 - unmarried_Ptilde_low)*Phi_approx(get_mu(
      unmarried_Ptilde_star,
      unmarried_Ptilde_raw_n_terms,
      n_geounit,
      unmarried_Ptilde_model_matrix_w, unmarried_Ptilde_model_matrix_v, unmarried_Ptilde_model_matrix_u));
  // omega
  vector[unmarried_Omega_raw_n_terms] unmarried_Omega_star = get_mu_star(
    unmarried_Omega_n_sigma, unmarried_Omega_n_sigma_fixed, unmarried_Omega_n_sigma_estimate,
    unmarried_Omega_sigma_fixed, unmarried_Omega_sigma_estimate,
    unmarried_Omega_scalarprior_mean, unmarried_Omega_scalarprior_sd,
    unmarried_Omega_raw_n_terms, unmarried_Omega_raw_n_terms_fixed, unmarried_Omega_raw_n_terms_estimate,
    unmarried_Omega_raw_fixed, unmarried_Omega_raw_estimate,
    unmarried_Omega_re_start, unmarried_Omega_re_end);
  vector[n_geounit] unmarried_Omega = get_mu(
    unmarried_Omega_star,
    unmarried_Omega_raw_n_terms,
    n_geounit,
    unmarried_Omega_model_matrix_w, unmarried_Omega_model_matrix_v, unmarried_Omega_model_matrix_u);

  // splines coeff
  matrix<lower=0>[unmarried_Betas_n_sigma_estimate, k] unmarried_Betas_sigma_estimate;
  if (unmarried_Betas_n_sigma_estimate > 0) {
    unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 1] = reverse(unmarried_Betas_sigma_estimate_reverse_1);
    unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 2] = reverse(unmarried_Betas_sigma_estimate_reverse_2);
    unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 3] = reverse(unmarried_Betas_sigma_estimate_reverse_3);
    unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 4] = reverse(unmarried_Betas_sigma_estimate_reverse_4);
    //unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 5] = reverse(unmarried_Betas_sigma_estimate_reverse_5);
    //   unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 6] = reverse(unmarried_Betas_sigma_estimate_reverse_6);
    //   unmarried_Betas_sigma_estimate[1:unmarried_Betas_n_sigma_estimate, 7] = reverse(unmarried_Betas_sigma_estimate_reverse_7);
  }
  matrix[unmarried_Betas_raw_n_terms,unmarried_Betas_k_terms] unmarried_Betas_star = get_mudimhk_star(k ,
       unmarried_Betas_n_sigma, unmarried_Betas_n_sigma_fixed, unmarried_Betas_n_sigma_estimate,
       unmarried_Betas_sigma_fixed, unmarried_Betas_sigma_estimate,
       unmarried_Betas_scalarprior_mean, unmarried_Betas_scalarprior_sd,
       unmarried_Betas_raw_n_terms, unmarried_Betas_raw_n_terms_fixed, unmarried_Betas_raw_n_terms_estimate,
       unmarried_Betas_raw_fixed, unmarried_Betas_raw_estimate,
       unmarried_Betas_re_start, unmarried_Betas_re_end);
   matrix[n_geounit,unmarried_Betas_k_terms] unmarried_Betas = get_mudimhk(k,
       unmarried_Betas_star,
       unmarried_Betas_raw_n_terms,
       n_geounit,
       unmarried_Betas_model_matrix_w, unmarried_Betas_model_matrix_v, unmarried_Betas_model_matrix_u);
  matrix[n_geounit,unmarried_Betas_k_terms] tr_unmarried_Betas_nonzero = unmarried_Betas_lower_bound + (unmarried_Betas_upper_bound - unmarried_Betas_lower_bound) * inv_logit(unmarried_Betas);


 // smoothing
  matrix[n_geounit, T] unmarried_Epsilon;
  unmarried_Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    unmarried_Rho_fixed, unmarried_Rho_estimate,
    unmarried_Tau_fixed, unmarried_Tau_estimate,
    unmarried_Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_unmarried_Eta_obs = process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing, unmarried_Epsilon,
                        unmarried_Omega, tr_unmarried_Ptilde, tr_unmarried_Betas_nonzero,
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

  // asymptote
  vector[unmarried_d_Ptilde_n_sigma_estimate] unmarried_d_Ptilde_sigma_estimate = reverse(unmarried_d_Ptilde_sigma_estimate_reverse);
  vector[unmarried_d_Ptilde_raw_n_terms] unmarried_d_Ptilde_star = get_mu_star(
      unmarried_d_Ptilde_n_sigma, unmarried_d_Ptilde_n_sigma_fixed, unmarried_d_Ptilde_n_sigma_estimate,
      unmarried_d_Ptilde_sigma_fixed, unmarried_d_Ptilde_sigma_estimate,
      unmarried_d_Ptilde_scalarprior_mean,
      unmarried_d_Ptilde_scalarprior_sd,
      unmarried_d_Ptilde_raw_n_terms, unmarried_d_Ptilde_raw_n_terms_fixed, unmarried_d_Ptilde_raw_n_terms_estimate,
      unmarried_d_Ptilde_raw_fixed, unmarried_d_Ptilde_raw_estimate,
      unmarried_d_Ptilde_re_start, unmarried_d_Ptilde_re_end);
//   vector[n_geounit] tr_unmarried_d_Ptilde = unmarried_d_Ptilde_low + (1 - unmarried_d_Ptilde_low)*inv_logit(get_mu(
// use probit
   vector[n_geounit] tr_unmarried_d_Ptilde = //rep_vector(0.9999, n_geounit);
   unmarried_d_Ptilde_low + (1 - unmarried_d_Ptilde_low)*Phi_approx(get_mu(
      unmarried_d_Ptilde_star,
      unmarried_d_Ptilde_raw_n_terms,
      n_geounit,
      unmarried_d_Ptilde_model_matrix_w, unmarried_d_Ptilde_model_matrix_v, unmarried_d_Ptilde_model_matrix_u));
  // omega
  vector[unmarried_d_Omega_raw_n_terms] unmarried_d_Omega_star = get_mu_star(
    unmarried_d_Omega_n_sigma, unmarried_d_Omega_n_sigma_fixed, unmarried_d_Omega_n_sigma_estimate,
    unmarried_d_Omega_sigma_fixed, unmarried_d_Omega_sigma_estimate,
    unmarried_d_Omega_scalarprior_mean, unmarried_d_Omega_scalarprior_sd,
    unmarried_d_Omega_raw_n_terms, unmarried_d_Omega_raw_n_terms_fixed, unmarried_d_Omega_raw_n_terms_estimate,
    unmarried_d_Omega_raw_fixed, unmarried_d_Omega_raw_estimate,
    unmarried_d_Omega_re_start, unmarried_d_Omega_re_end);
  vector[n_geounit] unmarried_d_Omega = get_mu(
    unmarried_d_Omega_star,
    unmarried_d_Omega_raw_n_terms,
    n_geounit,
    unmarried_d_Omega_model_matrix_w, unmarried_d_Omega_model_matrix_v, unmarried_d_Omega_model_matrix_u);

  // splines coeff
  matrix<lower=0>[unmarried_d_Betas_n_sigma_estimate, k] unmarried_d_Betas_sigma_estimate;
  if (unmarried_d_Betas_n_sigma_estimate > 0) {
    unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 1] = reverse(unmarried_d_Betas_sigma_estimate_reverse_1);
    unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 2] = reverse(unmarried_d_Betas_sigma_estimate_reverse_2);
    unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 3] = reverse(unmarried_d_Betas_sigma_estimate_reverse_3);
    unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 4] = reverse(unmarried_d_Betas_sigma_estimate_reverse_4);
    //unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 5] = reverse(unmarried_d_Betas_sigma_estimate_reverse_5);
    //   unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 6] = reverse(unmarried_d_Betas_sigma_estimate_reverse_6);
    //   unmarried_d_Betas_sigma_estimate[1:unmarried_d_Betas_n_sigma_estimate, 7] = reverse(unmarried_d_Betas_sigma_estimate_reverse_7);
  }
  matrix[unmarried_d_Betas_raw_n_terms,unmarried_d_Betas_k_terms] unmarried_d_Betas_star = get_mudimhk_star(k ,
       unmarried_d_Betas_n_sigma, unmarried_d_Betas_n_sigma_fixed, unmarried_d_Betas_n_sigma_estimate,
       unmarried_d_Betas_sigma_fixed, unmarried_d_Betas_sigma_estimate,
       unmarried_d_Betas_scalarprior_mean, unmarried_d_Betas_scalarprior_sd,
       unmarried_d_Betas_raw_n_terms, unmarried_d_Betas_raw_n_terms_fixed, unmarried_d_Betas_raw_n_terms_estimate,
       unmarried_d_Betas_raw_fixed, unmarried_d_Betas_raw_estimate,
       unmarried_d_Betas_re_start, unmarried_d_Betas_re_end);
   matrix[n_geounit,unmarried_d_Betas_k_terms] unmarried_d_Betas = get_mudimhk(k,
       unmarried_d_Betas_star,
       unmarried_d_Betas_raw_n_terms,
       n_geounit,
       unmarried_d_Betas_model_matrix_w, unmarried_d_Betas_model_matrix_v, unmarried_d_Betas_model_matrix_u);
  matrix[n_geounit,unmarried_d_Betas_k_terms] tr_unmarried_d_Betas_nonzero = unmarried_d_Betas_lower_bound + (unmarried_d_Betas_upper_bound - unmarried_d_Betas_lower_bound) * inv_logit(unmarried_d_Betas);


 // smoothing
  matrix[n_geounit, T] unmarried_d_Epsilon;
  unmarried_d_Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    unmarried_d_Rho_fixed, unmarried_d_Rho_estimate,
    unmarried_d_Tau_fixed, unmarried_d_Tau_estimate,
    unmarried_d_Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_unmarried_d_Eta_obs = process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing, unmarried_d_Epsilon,
                        unmarried_d_Omega, tr_unmarried_d_Ptilde, tr_unmarried_d_Betas_nonzero,
                        k, num_basis,   ext_knots,  spline_degree, add_shock);

  

  // z_Omega
  vector[z_Omega_raw_n_terms] z_Omega_star = get_mu_star(
    z_Omega_n_sigma, z_Omega_n_sigma_fixed, z_Omega_n_sigma_estimate,
    z_Omega_sigma_fixed, z_Omega_sigma_estimate,
    z_Omega_scalarprior_mean, z_Omega_scalarprior_sd,
    z_Omega_raw_n_terms, z_Omega_raw_n_terms_fixed, z_Omega_raw_n_terms_estimate,
    z_Omega_raw_fixed, z_Omega_raw_estimate,
    z_Omega_re_start, z_Omega_re_end);
  vector[n_geounit] z_Omega = get_mu(
    z_Omega_star,
    z_Omega_raw_n_terms,
    n_geounit,
    z_Omega_model_matrix_w, z_Omega_model_matrix_v, z_Omega_model_matrix_u);


 // smoothing
  matrix[n_geounit, T] z_Epsilon;
  z_Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, z_fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    z_Rho_fixed, z_Rho_estimate,
    z_Tau_fixed, z_Tau_estimate,
    z_Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_z_Eta_obs = rw1process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing,
                       z_Epsilon,
                        z_Omega, add_shock);



  // unmarried_z_Omega
  vector[unmarried_z_Omega_raw_n_terms] unmarried_z_Omega_star = get_mu_star(
    unmarried_z_Omega_n_sigma, unmarried_z_Omega_n_sigma_fixed, unmarried_z_Omega_n_sigma_estimate,
    unmarried_z_Omega_sigma_fixed, unmarried_z_Omega_sigma_estimate,
    unmarried_z_Omega_scalarprior_mean, unmarried_z_Omega_scalarprior_sd,
    unmarried_z_Omega_raw_n_terms, unmarried_z_Omega_raw_n_terms_fixed, unmarried_z_Omega_raw_n_terms_estimate,
    unmarried_z_Omega_raw_fixed, unmarried_z_Omega_raw_estimate,
    unmarried_z_Omega_re_start, unmarried_z_Omega_re_end);
  vector[n_geounit] unmarried_z_Omega = get_mu(
    unmarried_z_Omega_star,
    unmarried_z_Omega_raw_n_terms,
    n_geounit,
    unmarried_z_Omega_model_matrix_w, unmarried_z_Omega_model_matrix_v, unmarried_z_Omega_model_matrix_u);


 // smoothing
  matrix[n_geounit, T] unmarried_z_Epsilon;
  unmarried_z_Epsilon = compute_epsilon(
    smoothing, fix_subnat_corr, unmarried_z_fix_smoothing,
    rho_correlationeps_fixed, rho_correlationeps_estimate,
    unmarried_z_Rho_fixed, unmarried_z_Rho_estimate,
    unmarried_z_Tau_fixed, unmarried_z_Tau_estimate,
    unmarried_z_Epsilon_innovation, n_geounit, T, t_star,
    correlated_smoothing, max_cor_smoothing_block_size, n_cor_smoothing_blocks,
    cor_smoothing_block_sizes
  );

   matrix[n_geounit, T] tr_unmarried_z_Eta_obs = rw1process_model_returns_etatr(
                        n_geounit, T, t_star, t_min, t_max,
                        smoothing,
                       unmarried_z_Epsilon,
                        unmarried_z_Omega, add_shock);


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


  // minor note that this includes obs that are unmarried_NA, could subset instead
   array[add_dataoutliers * (1-unmarried_DM1_fix_nonse)] real<lower=0> unmarried_DM1_caux_dm_estimate = unmarried_DM1_sqrt_caux_dm_estimate^2;
   vector[unmarried_N] unmarried_DM1_scale;
   if (add_dataoutliers){
    unmarried_DM1_scale = get_scale(
        unmarried_DM1_local_shrinkage_dm, unmarried_DM1_nonse_fixed, unmarried_DM1_global_shrinkage_dm_fixed, unmarried_DM1_caux_dm_fixed,
        unmarried_DM1_nonse_estimate, unmarried_DM1_global_shrinkage_dm_estimate, unmarried_DM1_caux_dm_estimate,
        unmarried_N, unmarried_DM1_s, unmarried_isdhs, unmarried_source, unmarried_DM1_nooutlier, unmarried_S, unmarried_DM1_fix_nonse,
        unmarried_DM1_any_bias, unmarried_DM1_sdbias_fixed, unmarried_DM1_sdbias_estimate);
   } else {
    unmarried_DM1_scale = get_scale_nooutliers(
        unmarried_DM1_nonse_fixed,
        unmarried_DM1_nonse_estimate,
        unmarried_N, unmarried_DM1_s, unmarried_isdhs, unmarried_source, unmarried_DM1_nooutlier, unmarried_S, unmarried_DM1_fix_nonse);
   }
     // pma
  real<lower = 0, upper = 1> unmarried_DM1_rho_pma;
  if (unmarried_DM1_fix_nonse) {
    unmarried_DM1_rho_pma = unmarried_DM1_rho_pma_fixed[1];
  } else {
    unmarried_DM1_rho_pma = unmarried_DM1_rho_pma_estimate[1];
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


  // minor note that this includes obs that are unmarried_NA, could subset instead
   array[add_dataoutliers * (1-unmarried_DM2_fix_nonse)] real<lower=0> unmarried_DM2_caux_dm_estimate = unmarried_DM2_sqrt_caux_dm_estimate^2;
   vector[unmarried_N] unmarried_DM2_scale;
   if (add_dataoutliers){
    unmarried_DM2_scale = get_scale(
        unmarried_DM2_local_shrinkage_dm, unmarried_DM2_nonse_fixed, unmarried_DM2_global_shrinkage_dm_fixed, unmarried_DM2_caux_dm_fixed,
        unmarried_DM2_nonse_estimate, unmarried_DM2_global_shrinkage_dm_estimate, unmarried_DM2_caux_dm_estimate,
        unmarried_N, unmarried_DM2_s, unmarried_isdhs, unmarried_source, unmarried_DM2_nooutlier, unmarried_S, unmarried_DM2_fix_nonse,
        unmarried_DM2_any_bias, unmarried_DM2_sdbias_fixed, unmarried_DM2_sdbias_estimate);
   } else {
    unmarried_DM2_scale = get_scale_nooutliers(
        unmarried_DM2_nonse_fixed,
        unmarried_DM2_nonse_estimate,
        unmarried_N, unmarried_DM2_s, unmarried_isdhs, unmarried_source, unmarried_DM2_nooutlier, unmarried_S, unmarried_DM2_fix_nonse);
   }
     // pma
  real<lower = 0, upper = 1> unmarried_DM2_rho_pma;
  if (unmarried_DM2_fix_nonse) {
    unmarried_DM2_rho_pma = unmarried_DM2_rho_pma_fixed[1];
  } else {
    unmarried_DM2_rho_pma = unmarried_DM2_rho_pma_estimate[1];
  }
  // real<lower = 0, upper = 1> z_rho_pma;
  // if (add_z && fix_z) {
  //   z_rho_pma = z_rho_pma_fixed[1];
  // } else {
  //   z_rho_pma = z_rho_pma_estimate[1];
  // }


    // minor note that this includes obs that are NA, could subset instead
   array[add_dataoutliers * (1-DM3_fix_nonse)] real<lower=0> DM3_caux_dm_estimate = DM3_sqrt_caux_dm_estimate^2;
   vector[N] DM3_scale;
   if (add_dataoutliers){
    DM3_scale = get_scale(
        DM3_local_shrinkage_dm, DM3_nonse_fixed, DM3_global_shrinkage_dm_fixed, DM3_caux_dm_fixed,
        DM3_nonse_estimate, DM3_global_shrinkage_dm_estimate, DM3_caux_dm_estimate,
        N, DM3_s, isdhs, source, DM3_nooutlier, S, DM3_fix_nonse,
        DM3_any_bias, DM3_sdbias_fixed, DM3_sdbias_estimate);
   } else {
    DM3_scale = get_scale_nooutliers(
        DM3_nonse_fixed,
        DM3_nonse_estimate,
        N, DM3_s, isdhs, source, DM3_nooutlier, S, DM3_fix_nonse);
   }
     // pma
  real<lower = 0, upper = 1> DM3_rho_pma;
  if (DM3_fix_nonse) {
    DM3_rho_pma = DM3_rho_pma_fixed[1];
  } else {
    DM3_rho_pma = DM3_rho_pma_estimate[1];
  }
  // real<lower = 0, upper = 1> z_rho_pma;
  // if (add_z && fix_z) {
  //   z_rho_pma = z_rho_pma_fixed[1];
  // } else {
  //   z_rho_pma = z_rho_pma_estimate[1];
  // }


  // minor note that this includes obs that are unmarried_NA, could subset instead
   array[add_dataoutliers * (1-unmarried_DM3_fix_nonse)] real<lower=0> unmarried_DM3_caux_dm_estimate = unmarried_DM3_sqrt_caux_dm_estimate^2;
   vector[unmarried_N] unmarried_DM3_scale;
   if (add_dataoutliers){
    unmarried_DM3_scale = get_scale(
        unmarried_DM3_local_shrinkage_dm, unmarried_DM3_nonse_fixed, unmarried_DM3_global_shrinkage_dm_fixed, unmarried_DM3_caux_dm_fixed,
        unmarried_DM3_nonse_estimate, unmarried_DM3_global_shrinkage_dm_estimate, unmarried_DM3_caux_dm_estimate,
        unmarried_N, unmarried_DM3_s, unmarried_isdhs, unmarried_source, unmarried_DM3_nooutlier, unmarried_S, unmarried_DM3_fix_nonse,
        unmarried_DM3_any_bias, unmarried_DM3_sdbias_fixed, unmarried_DM3_sdbias_estimate);
   } else {
    unmarried_DM3_scale = get_scale_nooutliers(
        unmarried_DM3_nonse_fixed,
        unmarried_DM3_nonse_estimate,
        unmarried_N, unmarried_DM3_s, unmarried_isdhs, unmarried_source, unmarried_DM3_nooutlier, unmarried_S, unmarried_DM3_fix_nonse);
   }
     // pma
  real<lower = 0, upper = 1> unmarried_DM3_rho_pma;
  if (unmarried_DM3_fix_nonse) {
    unmarried_DM3_rho_pma = unmarried_DM3_rho_pma_fixed[1];
  } else {
    unmarried_DM3_rho_pma = unmarried_DM3_rho_pma_estimate[1];
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

  unmarried_DM1_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  unmarried_DM1_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    unmarried_DM1_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    unmarried_DM1_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    unmarried_DM1_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }

    DM2_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  DM2_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    DM2_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    DM2_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    DM2_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }

  unmarried_DM2_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  unmarried_DM2_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    unmarried_DM2_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    unmarried_DM2_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    unmarried_DM2_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }

    DM3_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  DM3_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    DM3_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    DM3_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    DM3_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
  }

  unmarried_DM3_sdbias_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];//2;
  unmarried_DM3_nonse_estimate ~ normal(0, 0.5)T[verysmallnumber, positive_infinity()];
  if (add_dataoutliers){
    unmarried_DM3_local_shrinkage_dm ~ student_t(1, 0, 1)T[0, positive_infinity()];
    unmarried_DM3_global_shrinkage_dm_estimate ~ student_t(1, 0, 0.04)T[0, positive_infinity()];
    unmarried_DM3_sqrt_caux_dm_estimate ~ normal(0,1)T[0, positive_infinity()];
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
  unmarried_Ptilde_raw_estimate ~ std_normal();
  unmarried_Omega_raw_estimate ~ std_normal();
  to_vector(unmarried_Betas_raw_estimate) ~ std_normal();

  // variances
  //unmarried_Ptilde_sigma_estimate ~ normal(0, unmarried_Ptilde_prior_sd_sigma_estimate)T[0, positive_infinity()];
  unmarried_Ptilde_sigma_estimate_reverse ~ normal(0, unmarried_Ptilde_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  unmarried_Omega_sigma_estimate ~ normal(0, unmarried_Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(unmarried_Betas_sigma_estimate) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[0, positive_infinity()];
   // one vector per spline coefficient
  to_vector(unmarried_Betas_sigma_estimate_reverse_1) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(unmarried_Betas_sigma_estimate_reverse_2) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(unmarried_Betas_sigma_estimate_reverse_3) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(unmarried_Betas_sigma_estimate_reverse_4) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(unmarried_Betas_sigma_estimate_reverse_5) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(unmarried_Betas_sigma_estimate_reverse_6) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(unmarried_Betas_sigma_estimate_reverse_7) ~ normal(0, unmarried_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(unmarried_Epsilon_innovation) ~ std_normal();
    if(fix_smoothing == 0) {
      unmarried_Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      unmarried_Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
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



  // hierarchical parameters
  unmarried_d_Ptilde_raw_estimate ~ std_normal();
  unmarried_d_Omega_raw_estimate ~ std_normal();
  to_vector(unmarried_d_Betas_raw_estimate) ~ std_normal();

  // variances
  //unmarried_d_Ptilde_sigma_estimate ~ normal(0, unmarried_d_Ptilde_prior_sd_sigma_estimate)T[0, positive_infinity()];
  unmarried_d_Ptilde_sigma_estimate_reverse ~ normal(0, unmarried_d_Ptilde_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  unmarried_d_Omega_sigma_estimate ~ normal(0, unmarried_d_Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(unmarried_d_Betas_sigma_estimate) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[0, positive_infinity()];
   // one vector per spline coefficient
  to_vector(unmarried_d_Betas_sigma_estimate_reverse_1) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(unmarried_d_Betas_sigma_estimate_reverse_2) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(unmarried_d_Betas_sigma_estimate_reverse_3) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  to_vector(unmarried_d_Betas_sigma_estimate_reverse_4) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  //to_vector(unmarried_d_Betas_sigma_estimate_reverse_5) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(unmarried_d_Betas_sigma_estimate_reverse_6) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];
  // to_vector(unmarried_d_Betas_sigma_estimate_reverse_7) ~ normal(0, unmarried_d_Betas_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(unmarried_d_Epsilon_innovation) ~ std_normal();
    if(fix_smoothing == 0) {
      unmarried_d_Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      unmarried_d_Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }



    // hierarchical parameters
  z_Omega_raw_estimate ~ std_normal();
  z_Omega_sigma_estimate ~ normal(0, z_Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(z_Epsilon_innovation) ~ std_normal();
    if(z_fix_smoothing == 0) {
      z_Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      z_Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }



  // hierarchical parameters
  unmarried_z_Omega_raw_estimate ~ std_normal();
  unmarried_z_Omega_sigma_estimate ~ normal(0, unmarried_z_Omega_prior_sd_sigma_estimate)T[verysmallnumber, positive_infinity()];

 // smoothing terms
  if(smoothing) {
    to_vector(unmarried_z_Epsilon_innovation) ~ std_normal();
    if(unmarried_z_fix_smoothing == 0) {
      unmarried_z_Rho_estimate[1] ~ normal(0, 0.5)T[0,1];
      unmarried_z_Tau_estimate[1] ~ normal(0, 0.2)T[verysmallnumber, positive_infinity()];
    }
  }




  // using one dm block here for now to hardcode what mean is
    // model
  // fit to the data

  // get means for all, to use for pma and non-pma
  vector[N] modern;
  vector[N] logit_unmetovernonmodern;
  vector[N] trad;
  vector[N] unmet;
  vector[N] logit_dm3;

  for (i in 1:N){
    if (geo_unit[i] == 0){ // aggregate obs, calculate national aggregates
      vector[n_geounit] demand_r = inv_tr_eta_colvector(tr_d_Eta_obs[, time[i]]);
      vector[n_geounit] modern_r = inv_tr_eta_colvector(tr_Eta_obs[, time[i]]) .* demand_r;
      vector[n_geounit] unmet_r = demand_r - modern_r;
      vector[n_geounit] trad_r = inv_tr_eta_colvector(tr_z_Eta_obs[, time[i]]) .* unmet_r;
      modern[i] = sum(modern_r .* to_vector(geo_unit_natpop_weight_tr[time[i], ]));
      unmet[i] = sum(unmet_r .* to_vector(geo_unit_natpop_weight_tr[time[i], ]));
      trad[i] = sum(trad_r .* to_vector(geo_unit_natpop_weight_tr[time[i], ]));
    } else {
      real demand = inv_tr_eta(tr_d_Eta_obs[geo_unit[i], time[i]]);
      modern[i] = inv_tr_eta(tr_Eta_obs[geo_unit[i], time[i]]) * demand;
      unmet[i] = demand - modern[i];
      // trad = z_eta .* unmet;
      trad[i] = inv_tr_eta(tr_z_Eta_obs[geo_unit[i], time[i]])*unmet[i];
    }
    logit_unmetovernonmodern[i] = logit(unmet[i]/(1-modern[i]));

    // for dm3, mean = ifelse(is.na(dm2), trad/unmet, trad/1minmodern)
    if (DM2_obs_isna[i] == 0){
      logit_dm3[i] = logit(trad[i]/unmet[i]);
    } else {
      logit_dm3[i] = logit(trad[i]/(1-modern[i]));
    }
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
      if(DM3_obs_isna[i] == 0) {
         DM3_y[i] ~ normal(logit_dm3[i], DM3_scale[i]);
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
    matrix[npma_c[index_country], npma_c[index_country]] DM3_cov_matrix_pma;
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
        DM3_cov_matrix_pma[j1, j2] =
          DM3_scale[j1_index]* DM3_scale[j2_index] *
            DM3_rho_pma^(abs(exact_time[j1_index] - exact_time[j2_index]));

    }}
    DM1_y[indices_inc] ~ multi_normal(logit(modern[indices_inc]), DM1_cov_matrix_pma);
    DM2_y[indices_inc] ~ multi_normal(logit_unmetovernonmodern[indices_inc], DM2_cov_matrix_pma);
    DM3_y[indices_inc] ~ multi_normal(logit_dm3[indices_inc], DM3_cov_matrix_pma);
  }



  // get means for all, to use for pma and non-pma
   // get means for all, to use for pma and non-pma
  vector[unmarried_N] unmarried_modern;
  vector[unmarried_N] unmarried_logit_unmetovernonmodern;
  vector[unmarried_N] unmarried_trad;
  vector[unmarried_N] unmarried_unmet;
  vector[unmarried_N] unmarried_logit_dm3;

  for (i in 1:unmarried_N){
    if (unmarried_geo_unit[i] == 0){ // aggregate obs, calculate national aggregates
      vector[n_geounit] unmarried_demand_r = inv_tr_eta_colvector(tr_unmarried_d_Eta_obs[, unmarried_time[i]]);
      vector[n_geounit] unmarried_modern_r = inv_tr_eta_colvector(tr_unmarried_Eta_obs[, unmarried_time[i]]) .* unmarried_demand_r;
      vector[n_geounit] unmarried_unmet_r = unmarried_demand_r - unmarried_modern_r;
      vector[n_geounit] unmarried_trad_r = inv_tr_eta_colvector(tr_unmarried_z_Eta_obs[, unmarried_time[i]]) .* unmarried_unmet_r;
      unmarried_modern[i] = sum(unmarried_modern_r .* to_vector(unmarried_geo_unit_natpop_weight_tr[unmarried_time[i], ]));
      unmarried_unmet[i] = sum(unmarried_unmet_r .* to_vector(unmarried_geo_unit_natpop_weight_tr[unmarried_time[i], ]));
      unmarried_trad[i] = sum(unmarried_trad_r .* to_vector(unmarried_geo_unit_natpop_weight_tr[unmarried_time[i], ]));
    } else {
      real unmarried_demand = inv_tr_eta(tr_unmarried_d_Eta_obs[unmarried_geo_unit[i], unmarried_time[i]]);
      unmarried_modern[i] = inv_tr_eta(tr_unmarried_Eta_obs[unmarried_geo_unit[i], unmarried_time[i]]) * unmarried_demand;
      unmarried_unmet[i] = unmarried_demand - unmarried_modern[i];
      // trad = z_eta .* unmet;
      unmarried_trad[i] = inv_tr_eta(tr_unmarried_z_Eta_obs[unmarried_geo_unit[i], unmarried_time[i]]) * unmarried_unmet[i];
    }
    unmarried_logit_unmetovernonmodern[i] = logit(unmarried_unmet[i]/(1-unmarried_modern[i]));

    // for dm3, mean = ifelse(is.na(dm2), trad/unmet, trad/1minmodern)
    if (unmarried_DM2_obs_isna[i] == 0){
      unmarried_logit_dm3[i] = logit(unmarried_trad[i]/unmarried_unmet[i]);
    } else {
      unmarried_logit_dm3[i] = logit(unmarried_trad[i]/(1-unmarried_modern[i]));
    }
  }

  // non-pma data
  for(i in 1:unmarried_N) {
    if(unmarried_held_out[i] == 0 && unmarried_ispma[i] == 0) {
      if(unmarried_DM1_obs_isna[i] == 0) {
          unmarried_DM1_y[i] ~ normal(logit(unmarried_modern[i]), unmarried_DM1_scale[i]);
      }
      if(unmarried_DM2_obs_isna[i] == 0) {
         unmarried_DM2_y[i] ~ normal(unmarried_logit_unmetovernonmodern[i], unmarried_DM2_scale[i]);
      }
      if(unmarried_DM3_obs_isna[i] == 0) {
         unmarried_DM3_y[i] ~ normal(unmarried_logit_dm3[i], unmarried_DM3_scale[i]);
      }
    }
  }
  // pma
  // no densities, constrained between 0 and 1
  // not yet for held out (held out just filtered out)
  // not yet for agregates
  // same for unmet
  for (index_country in unmarried_pma_country_indices){
    array[unmarried_npma_c[index_country]] int unmarried_indices_inc =
                unmarried_indices_pma[unmarried_nstart_pma_c[index_country]:
          (unmarried_nstart_pma_c[index_country] + unmarried_npma_c[index_country] - 1)];
    matrix[unmarried_npma_c[index_country], unmarried_npma_c[index_country]] unmarried_DM1_cov_matrix_pma;
    matrix[unmarried_npma_c[index_country], unmarried_npma_c[index_country]] unmarried_DM2_cov_matrix_pma;
    matrix[unmarried_npma_c[index_country], unmarried_npma_c[index_country]] unmarried_DM3_cov_matrix_pma;
    for (j1 in 1:unmarried_npma_c[index_country]){
      for (j2 in 1:unmarried_npma_c[index_country]){
        int unmarried_j1_index = unmarried_indices_inc[j1];
        int unmarried_j2_index = unmarried_indices_inc[j2];
        unmarried_DM1_cov_matrix_pma[j1, j2] =
            unmarried_DM1_scale[unmarried_j1_index]* unmarried_DM1_scale[unmarried_j2_index] *
              unmarried_DM1_rho_pma^(abs(unmarried_exact_time[unmarried_j1_index] -
              unmarried_exact_time[unmarried_j2_index]));
        unmarried_DM2_cov_matrix_pma[j1, j2] =
            unmarried_DM2_scale[unmarried_j1_index]* unmarried_DM2_scale[unmarried_j2_index] *
              unmarried_DM2_rho_pma^(abs(unmarried_exact_time[unmarried_j1_index] -
              unmarried_exact_time[unmarried_j2_index]));
        unmarried_DM3_cov_matrix_pma[j1, j2] =
          unmarried_DM3_scale[unmarried_j1_index]* unmarried_DM3_scale[unmarried_j2_index] *
            unmarried_DM3_rho_pma^(abs(unmarried_exact_time[unmarried_j1_index] -
            unmarried_exact_time[unmarried_j2_index]));

    }}
    unmarried_DM1_y[unmarried_indices_inc] ~ multi_normal(logit(unmarried_modern[unmarried_indices_inc]),
                      unmarried_DM1_cov_matrix_pma);
    unmarried_DM2_y[unmarried_indices_inc] ~ multi_normal(unmarried_logit_unmetovernonmodern[unmarried_indices_inc],
                      unmarried_DM2_cov_matrix_pma);
    unmarried_DM3_y[unmarried_indices_inc] ~ multi_normal(unmarried_logit_dm3[unmarried_indices_inc],
                      unmarried_DM3_cov_matrix_pma);
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

  matrix[n_geounit, T] unmarried_Eta = inv_tr_eta_matrix(
            process_model_outsideobs(tr_unmarried_Eta_obs,
                      n_geounit, T,  t_min, t_max,
                      smoothing, unmarried_Epsilon,
                      tr_unmarried_Ptilde, tr_unmarried_Betas_nonzero,
                      k, num_basis,   ext_knots,  spline_degree, add_shock));

    matrix[n_geounit, T] d_Eta = inv_tr_eta_matrix(
            process_model_outsideobs(tr_d_Eta_obs,
                      n_geounit, T,  t_min, t_max,
                      smoothing, d_Epsilon,
                      tr_d_Ptilde, tr_d_Betas_nonzero,
                      k, num_basis,   ext_knots,  spline_degree, add_shock));

  matrix[n_geounit, T] unmarried_d_Eta = inv_tr_eta_matrix(
            process_model_outsideobs(tr_unmarried_d_Eta_obs,
                      n_geounit, T,  t_min, t_max,
                      smoothing, unmarried_d_Epsilon,
                      tr_unmarried_d_Ptilde, tr_unmarried_d_Betas_nonzero,
                      k, num_basis,   ext_knots,  spline_degree, add_shock));

  
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

  // allwomen results
// we have unmet, trad, and d_Eta (demand), for "" and "unmarried_"
// results per geounit
matrix[n_geounit, T] allwomen_unmet = prop_married_rt .* unmet +
                                        (1-prop_married_rt) .* unmarried_unmet;
matrix[n_geounit, T] allwomen_trad = prop_married_rt .* trad +
                                        (1-prop_married_rt) .* unmarried_trad;
matrix[n_geounit, T] allwomen_demand = prop_married_rt .* d_Eta +
                                        (1-prop_married_rt) .* unmarried_d_Eta;
matrix[n_geounit, T] allwomen_modern = allwomen_demand - allwomen_unmet;




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



  // subnational aggregates
  // assume we only do this with unmarried_traditional, then we have unmarried_trad and unmarried_unmet
  // (and unmarried_unmet, so could update this)
  // vector[T] unmarried_demand_aggr;
  vector[T] unmarried_modern_aggr;
  vector[T] unmarried_unmet_aggr;
  vector[T] unmarried_trad_aggr;
  for(t in 1:T) {
    // weighting is for the proportions (so unmarried_demand, need to get unmarried_modern first at regional level)
     // unmarried_demand_aggr[t] = sum(unmarried_d_Eta[,t] .* to_vector(unmarried_geo_unit_natpop_weight_tr[t, ]));
     //vector[n_geounit] unmarried_modern_regions = unmarried_d_Eta[,t] .* Eta[,t];
     unmarried_modern_aggr[t] = sum(unmarried_modern[,t] .* to_vector(unmarried_geo_unit_natpop_weight_tr[t, ]));
     unmarried_unmet_aggr[t] = sum(unmarried_unmet[, t] .* to_vector(unmarried_geo_unit_natpop_weight_tr[t, ]));
     unmarried_trad_aggr[t] = sum(unmarried_trad[, t] .* to_vector(unmarried_geo_unit_natpop_weight_tr[t, ]));
  }



  
// also add national level ALL women results
// assuming we don't run aggregates without running all-women
// we have allwomen_unmet, trad, demand, modern per geo_unit
// now get national level

vector[T] allwomen_modern_aggr= prop_married_t .* modern_aggr +
                                        (1-prop_married_t) .* unmarried_modern_aggr;
vector[T] allwomen_unmet_aggr= prop_married_t .* unmet_aggr +
                                        (1-prop_married_t) .* unmarried_unmet_aggr;
vector[T] allwomen_trad_aggr= prop_married_t .* trad_aggr +
                                        (1-prop_married_t) .* unmarried_trad_aggr;



}






