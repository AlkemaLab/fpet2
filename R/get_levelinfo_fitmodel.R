
#' Get list with information on hierarchical levels for transition model parameters,
#' and which ones to fix
#'
#' @param runstep character, from model fit
#' @param hierarchical_terms list of hierarchical terms from fit_model inputs or NULL
#' @param global_fit fit object from a global fit
#' @param prefix character, prefix for process model, e.g., "", d_", "z_"
#' @param add_subnational_hierarchy character, name of subnational hierarchy level to add, e.g., "subnat"
#'
#' @returns list with hierarchical terms and what's fixed
#' @keywords internal
get_hierlevels <- function(runstep,
                           hierarchical_terms = NULL,
                           global_fit = NULL,
                           prefix = "",
                           add_subnational_hierarchy = NULL){
  # returns hier levels and _fixed versions
  # no prefixes needed here, as long as it's saved in a list with name that indicates the process model

  if (runstep %in% c("step1a")){
    # take levels from transitionmodelparam
    print("For 1a, we take all levels from the fit_model inputs")
    hierarchical_terms_and_fixed = hierarchical_terms
    # and don't fix all
    print("We do not fix any terms or sigmas of hierarchical models.")
    # actually these fixed terms not used in localhierarchy::hierarchical_param_stan_data
    # hierarchical_terms_and_fixed[["hierarchical_level_sigmas_fixed"]] = c()
    # hierarchical_terms_and_fixed$hierarchical_level_terms_fixed = c()
    # hierarchical_terms_and_fixed$hierarchical_asymptote_sigmas_fixed = c()
    # hierarchical_terms_and_fixed$hierarchical_asymptote_terms_fixed = c()
    # hierarchical_terms_and_fixed$hierarchical_splines_sigmas_fixed = c()
    # hierarchical_terms_and_fixed$hierarchical_splines_terms_fixed = c()
    return(hierarchical_terms_and_fixed)
  }
  # get this outside list first, then create list
  print("We take all hier terms from the global fit, using prefix")
  hierarchical_level = global_fit[[paste0(prefix, "hierarchical_terms_and_fixed")]]$hierarchical_level
  hierarchical_splines = global_fit[[paste0(prefix, "hierarchical_terms_and_fixed")]]$hierarchical_splines
  hierarchical_asymptote = global_fit[[paste0(prefix, "hierarchical_terms_and_fixed")]]$hierarchical_asymptote
  # hierarchical_level = global_fit$hierarchical_level
  # hierarchical_splines = global_fit$hierarchical_splines
  # hierarchical_asymptote = global_fit$hierarchical_asymptote
  # update june 24: add local_subnational, extra level needs to be added
  if (runstep %in% c("global_subnational", "local_subnational")){
    # note: for local_subnational, we don't add a level!
    print("For subnational global run, we add a level for subnational hierarchical settings ")
    hierarchical_splines <- c(hierarchical_splines, add_subnational_hierarchy)
    hierarchical_level <- c(hierarchical_level, add_subnational_hierarchy)
    hierarchical_asymptote <- c(hierarchical_asymptote, add_subnational_hierarchy)
  }
  # consider what to fix in hierarchical set up
  # for mean terms, we always fit up to the 2nd lowest level
  print("For hierarchical terms, we fix things up to the 2nd-lowest level.")
  # could consider for subnational global, to not fit country means
  hierarchical_asymptote_terms_fixed = hierarchical_asymptote[1:(length(hierarchical_asymptote)-1)]
  hierarchical_splines_terms_fixed = hierarchical_splines[1:(length(hierarchical_splines)-1)]
  hierarchical_level_terms_fixed = hierarchical_level[1:(length(hierarchical_level)-1)]
  # for sigma, differs between local and not local, and step 2 versus the other steps!
  if (runstep %in% c("local_national", "local_subnational", "step2")){
    print("We fix all sigmas of hierarchical models for demand and ds.")
    hierarchical_asymptote_sigmas_fixed = hierarchical_asymptote[1:(length(hierarchical_asymptote))]
    hierarchical_splines_sigmas_fixed = hierarchical_splines[1:(length(hierarchical_splines))]
    hierarchical_level_sigmas_fixed = hierarchical_level[1:(length(hierarchical_level))]
  } else {
    print("For sigma terms in hierarchical models for demand and ds, we fix things up to the 2nd-lowest level.")
    # this includes the intercept!
    hierarchical_asymptote_sigmas_fixed = hierarchical_asymptote[1:(length(hierarchical_asymptote)-1)]
    hierarchical_splines_sigmas_fixed = hierarchical_splines[1:(length(hierarchical_splines)-1)]
    hierarchical_level_sigmas_fixed = hierarchical_level[1:(length(hierarchical_level)-1)]
  }

  hierarchical_terms_and_fixed = list(
    hierarchical_level  =  hierarchical_level,
    hierarchical_splines  = hierarchical_splines,
    hierarchical_asymptote = hierarchical_asymptote,
    hierarchical_level_sigmas_fixed = hierarchical_level_sigmas_fixed,
    hierarchical_level_terms_fixed = hierarchical_level_terms_fixed,
    hierarchical_asymptote_sigmas_fixed = hierarchical_asymptote_sigmas_fixed,
    hierarchical_asymptote_terms_fixed = hierarchical_asymptote_terms_fixed,
    hierarchical_splines_sigmas_fixed = hierarchical_splines_sigmas_fixed,
    hierarchical_splines_terms_fixed = hierarchical_splines_terms_fixed)

  # some checks re what's fixed:
  # It is not valid to fix terms at a given hierarchy level without also fixing
  # the sigma estimate at that hierarchy level.
  # For example, if x = mu + z * sigma,
  # it does not make sense to fix z without also fixing sigma.
  # on the other hand, we might fix sigma but not z if we want to borrow information
  # about variability from global fit, but the global fit didn't produce an estimate
  # of z for the geo unit we're interested in, or we are ok with re-estimating it?
  if (!all(hierarchical_asymptote_terms_fixed %in% hierarchical_asymptote_sigmas_fixed)) {
    stop("All values of hierarchical_asymptote_terms_fixed must also be contained in hierarchical_asymptote_sigmas_fixed")
  }
  if (!all(hierarchical_level_terms_fixed %in% hierarchical_level_sigmas_fixed)) {
    stop("All values of hierarchical_level_terms_fixed must also be contained in hierarchical_level_sigmas_fixed")
  }
  if (!all(hierarchical_splines_terms_fixed %in% hierarchical_splines_sigmas_fixed)) {
    stop("All values of hierarchical_splines_terms_fixed must also be contained in hierarchical_splines_sigmas_fixed")
  }

  return(hierarchical_terms_and_fixed)
}


#' Get list with information on hierarchical levels and which ones to fix for omegatrad
#'
#' @param runstep character, from model fit
#' @param hierarchical_terms list of hierarchical terms from fit_model inputs or NULL
#' @param global_fit fit object from a global fit
#' @param prefix  character, prefix for process model (should be "z_" for traditional)
#' @param add_subnational_hierarchy character, name of subnational hierarchy level to add, e.g., "subnat"
#'
#' @returns list with hierarchical terms and what's fixed
#' @keywords internal
get_hierlevels_omegatrad <- function(runstep,
                           hierarchical_terms = NULL,
                           global_fit = NULL,
                           prefix = "",
                           add_subnational_hierarchy = NULL){
  # returns hier levels and _fixed versions
  # no prefixes needed here, as long as it's saved in a list with name that indicates the process model

  if (runstep %in% c("step2")){
    # take levels from transitionmodelparam
    print("For 2, we take all levels for trad from the fit_model inputs")
    hierarchical_terms_and_fixed = hierarchical_terms
    # and don't fix all
    print("We do not fix any terms or sigmas of hierarchical models.")
    return(hierarchical_terms_and_fixed)
  }
  # get this outside list first, then create list
  print("We take all hier terms from the global fit, using prefix")
  hierarchical_level = global_fit[[paste0(prefix, "hierarchical_terms_and_fixed")]]$hierarchical_level
  # if (runstep %in% c("global_subnational")){
  #   # note: for local_subnational, we don't add a level!
  #   print("For subnational global run, we add a level for subnational hierarchical settings ")
  #   hierarchical_level <- c(hierarchical_level, add_subnational_hierarchy)
  # }
  # update june 24: add local_subnational, extra level needs to be added
  if (runstep %in% c("global_subnational", "local_subnational")){
    print("For subnational global run, and subnational from national, we add a level for subnational hierarchical settings ")
    hierarchical_level <- c(hierarchical_level, add_subnational_hierarchy)
  }
  # consider what to fix in hierarchical set up
  # for mean terms, we always fit up to the 2nd lowest level
  print("For hierarchical terms, we fix things up to the 2nd-lowest level.")
  hierarchical_level_terms_fixed = hierarchical_level[1:(length(hierarchical_level)-1)]
  # for sigma, differs between local and not local
  if (runstep %in% c("local_national", "local_subnational")){
    print("We fix all sigmas of hierarchical models.")
    hierarchical_level_sigmas_fixed = hierarchical_level[1:(length(hierarchical_level))]
  } else {
    print("For sigma terms in hierarchical models, we fix things up to the 2nd-lowest level.")
    hierarchical_level_sigmas_fixed = hierarchical_level[1:(length(hierarchical_level)-1)]
  }

  hierarchical_terms_and_fixed = list(
    hierarchical_level  =  hierarchical_level,
    hierarchical_level_sigmas_fixed = hierarchical_level_sigmas_fixed,
    hierarchical_level_terms_fixed = hierarchical_level_terms_fixed)

  # some checks re what's fixed:
  # It is not valid to fix terms at a given hierarchy level without also fixing
  # the sigma estimate at that hierarchy level.
  # For example, if x = mu + z * sigma,
  # it does not make sense to fix z without also fixing sigma.
  # on the other hand, we might fix sigma but not z if we want to borrow information
  # about variability from global fit, but the global fit didn't produce an estimate
  # of z for the geo unit we're interested in, or we are ok with re-estimating it?

  if (!all(hierarchical_level_terms_fixed %in% hierarchical_level_sigmas_fixed)) {
    stop("All values of hierarchical_level_terms_fixed must also be contained in hierarchical_level_sigmas_fixed")
  }


  return(hierarchical_terms_and_fixed)
}

