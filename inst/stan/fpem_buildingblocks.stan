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
  // int<lower=0, upper=1> DM1_fix_nonse;
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

  {{DATA_DATA}}

  // for parameters of process model
  // in old code, "d_" referred to demand
  // here do ds
  //modern[m] = eta[m] .* d_eta[m];
  //unmet[m] = d_eta[m] - modern[m];

  {{PROCESSMODEL_DATA}}
  {{PROCESSMODEL_DEMAND_DATA}}
  {{PROCESSMODEL_TRAD_DATA}}



  {{DATAMODEL_DATA}}
  {{DATAMODEL2_DATA}}
  {{DATAMODEL3_DATA}}

  {{AGGREGATES_DATA}}

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

     {{EMU_DATA}}

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
  {{PROCESSMODEL_TRANSFORMEDDATA}}
  {{PROCESSMODEL_DEMAND_TRANSFORMEDDATA}}
  {{PROCESSMODEL_TRAD_TRANSFORMEDDATA}}


}

/////////////////////////////////////////////////////
parameters {

  // for parameters of process model
  {{PROCESSMODEL_PARAMETERS}}
  {{PROCESSMODEL_DEMAND_PARAMETERS}}
   {{PROCESSMODEL_TRAD_PARAMETERS}}

  // Data model, one set per dm
  {{DATAMODEL_PARAMETERS}}
  {{DATAMODEL2_PARAMETERS}}
   {{DATAMODEL3_PARAMETERS}}

  //{{subnational_parameters}}
  // for subnational
  array[fix_subnat_corr ? 0 : 1] real<lower=0, upper = 1> rho_correlationeps_estimate;   // for correlated eps
   {{EMU_PARAMETERS}}

}

/////////////////////////////////////////////////////
transformed parameters {

  {{PROCESSMODEL_TRANSFORMEDPARAMETERS}}
  {{PROCESSMODEL_DEMAND_TRANSFORMEDPARAMETERS}}
  {{PROCESSMODEL_TRAD_TRANSFORMEDPARAMETERS}}

  {{DATAMODEL_TRANSFORMEDPARAMETERS}}
  {{DATAMODEL2_TRANSFORMEDPARAMETERS}}
  {{DATAMODEL3_TRANSFORMEDPARAMETERS}}

}

/////////////////////////////////////////////////////
model {

  // data model parameters, one set per dm
  {{DATAMODEL_MODEL_PARAM}}
  {{DATAMODEL2_MODEL_PARAM}}
  {{DATAMODEL3_MODEL_PARAM}}

  {{PROCESSMODEL_MODEL}}
  {{PROCESSMODEL_DEMAND_MODEL}}
  {{PROCESSMODEL_TRAD_MODEL}}

  // using one dm block here for now to hardcode what mean is
  {{DATAMODEL_MODEL}}

  // subnational
   if (!fix_subnat_corr){
    rho_correlationeps_estimate[1] ~ uniform(0, 1);
    }

    {{EMU_MODEL}}

}


generated quantities {

  {{PROCESSMODEL_GENERATEDQUANTITIES}}
  {{PROCESSMODEL_DEMAND_GENERATEDQUANTITIES}}
  {{PROCESSMODEL_TRAD_GENERATEDQUANTITIES}}
  {{ALLWOMEN_GENERATEDQUANTITIES}}
  {{AGGREGATES_GENERATEDQUANTITIES}}
  {{AGGREGATESALLWOMEN_GENERATEDQUANTITIES}}

}





