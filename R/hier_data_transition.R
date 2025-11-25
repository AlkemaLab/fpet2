
#' Get hierarchical data for transition model
#'
#' @param geo_unit_index geo_unit_index for model fit
#' @param prefix prefix to add to parameter names, e.g., "d_" or "z_"
#' @param stan_data_settings list of stan data settings, e.g., bounds and priors
#' @param stan_spline_data stan data for splines, needed for dimension of splines coeff vector
#' @param hierarchical_terms_and_fixed list with info on levels and what's fixed from other function
#' @param global_fit fit object from a global fit, or NULL if not doing a local fit
#'
#' @returns list with hier_stan_data_settings, hier_stan_data, and hier_data
#' @keywords internal
#'
get_hier_data_transition <- function(geo_unit_index,
                                     prefix = "",
                                     stan_data_settings,
                                     stan_spline_data, # needed for dimension of splines coeff vector
                                    hierarchical_terms_and_fixed, #list with info
                                    # info on levels and what's fixed from other function
                                    global_fit = NULL
                                     ){


  # stan_data_settings = list(
  #   Betas_upper_bound = 0.5,
  #   Betas_lower_bound = 0.01,
  #   Ptilde_low = 0,
  #   Ptilde_scalarprior_sd = 2,
  #   Ptilde_scalarprior_mean = 2,
  #   Ptilde_prior_sd_sigma_estimate = 1,
  #   Omega_scalarprior_sd = 2,
  #   Omega_scalarprior_mean = 0,
  #   Omega_prior_sd_sigma_estimate = 1,
  #   Betas_scalarprior_sd = 2,
  #   Betas_scalarprior_mean = -1,
  #   Betas_prior_sd_sigma_estimate = 1))

  hier_data <- list()
  if (is.null(global_fit)){
    # add prefix name
    #tmp <- list(x = 1, y = 2)
    #names(tmp) <- paste0("d_", names(tmp))
    #tmp
    names(stan_data_settings) <- paste0(prefix, names(stan_data_settings))
    hier_stan_data_settings <- stan_data_settings
  } else {
    # we do need to add the prefix in the final thing...
    hier_stan_data_settings <- global_fit[[paste0(prefix, "hier_stan_data_settings")]]
  }

  hier_stan_data <- hier_stan_data_settings
  hier_stan_data[[paste0(prefix, "Ptilde_isvector")]] <- FALSE
  hier_stan_data[[paste0(prefix, "Omega_isvector")]] <- FALSE
  hier_stan_data[[paste0(prefix, "Betas_isvector")]] <- TRUE

  parname <- paste0(prefix, "Ptilde")
  hier_data[[paste0(parname, "_data")]] <- localhierarchy::hierarchical_data(geo_unit_index,
                                                  hierarchical_terms_and_fixed$hierarchical_asymptote)
  hier_stan_data[[parname]] <- localhierarchy::hierarchical_param_stan_data(
    global_fit = global_fit,
    param_name = parname,
    param_data = hier_data[[paste0(parname, "_data")]],
    hierarchical_terms_fixed = hierarchical_terms_and_fixed$hierarchical_asymptote_terms_fixed,
    hierarchical_sigmas_fixed = hierarchical_terms_and_fixed$hierarchical_asymptote_sigmas_fixed)
  parname <- paste0(prefix, "Omega")
  hier_data[[paste0(parname, "_data")]] <- localhierarchy::hierarchical_data(geo_unit_index,
                                                 hierarchical_terms_and_fixed$hierarchical_level)
  hier_stan_data[[parname]] <- localhierarchy::hierarchical_param_stan_data(
    global_fit = global_fit,
    param_name = parname,
    param_data = hier_data[[paste0(parname, "_data")]],
    hierarchical_terms_fixed = hierarchical_terms_and_fixed$hierarchical_level_terms_fixed,
    hierarchical_sigmas_fixed = hierarchical_terms_and_fixed$hierarchical_level_sigmas_fixed)

  parname <- paste0(prefix, "Betas")
  # k is being calculated somewhere.... replaces this
  hier_stan_data[[paste0(parname, "_k_terms")]] <- stan_spline_data[["k"]]
  hier_data[[paste0(parname, "_data")]] <- localhierarchy::hierarchical_data(geo_unit_index,
                                                 hierarchical_terms_and_fixed$hierarchical_splines)
  hier_stan_data[[parname]] <- localhierarchy::hierarchical_param_stan_data(
    global_fit = global_fit,
    param_name = parname,
    param_data = hier_data[[paste0(parname, "_data")]],
    hierarchical_terms_fixed = hierarchical_terms_and_fixed$hierarchical_splines_terms_fixed,
    hierarchical_sigmas_fixed = hierarchical_terms_and_fixed$hierarchical_splines_sigmas_fixed)

  hier_stan_data <- purrr::list_flatten(hier_stan_data, name_spec = "{inner}")

  return(list(hier_stan_data_settings = hier_stan_data_settings,
              hier_stan_data = hier_stan_data,
              hier_data = hier_data))

  #return(append(hier_stan_data_settings,
  #        hier_stan_data,
  #        hier_data))

}



#' Get hierarchical data for Omegas in transition model
#'
#' Perhaps to combine with `get_hier_data_transition`
#'
#' @param geo_unit_index geo_unit_index for model fit
#' @param prefix prefix to add to parameter names, e.g., "d_" or "z_"
#' @param stan_data_settings list of stan data settings, e.g., bounds and priors
#' @param hierarchical_terms_and_fixed list with info on levels and what's fixed from other function
#' @param global_fit fit object from a global fit, or NULL if not doing a local fit
#' @param runstep what model is fitted
#'
#' @returns list with hier_stan_data_settings, hier_stan_data, and hier_data
#' @keywords internal
get_hier_data_transition_justomega <- function(geo_unit_index,
                                     prefix = "",
                                     stan_data_settings,
                                     hierarchical_terms_and_fixed, #list with info
                                     # info on levels and what's fixed from other function
                                     global_fit = NULL,
                                     runstep
){
  hier_data <- list()
  if (runstep == "step2"){
    names(stan_data_settings) <- paste0(prefix, names(stan_data_settings))
    hier_stan_data_settings <- stan_data_settings
  } else {
    # we do need to add the prefix in the final thing...
    hier_stan_data_settings <- global_fit[[paste0(prefix, "hier_stan_data_settings")]]
  }
  hier_stan_data <- hier_stan_data_settings
  hier_stan_data[[paste0(prefix, "Omega_isvector")]] <- FALSE
  parname <- paste0(prefix, "Omega")
  hier_data[[paste0(parname, "_data")]] <- localhierarchy::hierarchical_data(geo_unit_index,
                                                             hierarchical_terms_and_fixed$hierarchical_level)
  # avoid reliance on global fit non null argument...
  if (runstep == "step2"){
    hier_stan_data[[parname]] <- localhierarchy::hierarchical_param_stan_data(
      global_fit = NULL,
      param_name = parname,
      param_data = hier_data[[paste0(parname, "_data")]],
      hierarchical_terms_fixed = hierarchical_terms_and_fixed$hierarchical_level_terms_fixed,
      hierarchical_sigmas_fixed = hierarchical_terms_and_fixed$hierarchical_level_sigmas_fixed)
  } else {
    hier_stan_data[[parname]] <- localhierarchy::hierarchical_param_stan_data(
      global_fit = global_fit,
      param_name = parname,
      param_data = hier_data[[paste0(parname, "_data")]],
      hierarchical_terms_fixed = hierarchical_terms_and_fixed$hierarchical_level_terms_fixed,
      hierarchical_sigmas_fixed = hierarchical_terms_and_fixed$hierarchical_level_sigmas_fixed)
  }

  hier_stan_data <- purrr::list_flatten(hier_stan_data, name_spec = "{inner}")

  return(list(hier_stan_data_settings = hier_stan_data_settings,
              hier_stan_data = hier_stan_data,
              hier_data = hier_data))

}
# not used
# # if there are maxes for a sigma, update to get max for that sigma
# if (!is.null(global_fit)){
#   hier_stan_data_sigmamax <- list(
#     Ptilde_sigma_max = min(hier_stan_data$Ptilde_sigma_fixed),
#     Omega_sigma_max = min(hier_stan_data$Omega_sigma_fixed),
#     Betas_sigma_max_1 = min(hier_stan_data$Betas_sigma_fixed[,1]),
#     Betas_sigma_max_2 = min(hier_stan_data$Betas_sigma_fixed[,2]),
#     Betas_sigma_max_3 = min(hier_stan_data$Betas_sigma_fixed[,3]),
#     Betas_sigma_max_4 = min(hier_stan_data$Betas_sigma_fixed[,4])
#   )
# } else {
#   hier_stan_data_sigmamax <- list(
#     Ptilde_sigma_max = 5,
#     Omega_sigma_max = 5,
#     Betas_sigma_max_1  = 5,
#     Betas_sigma_max_2 = 5,
#     Betas_sigma_max_3 = 5,
#     Betas_sigma_max_4 = 5
#   )
# }
