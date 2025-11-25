
#' Plot hierarchical checks
#'
#' @param fit  fit object from a global fit
#' @param is_married logical, is this a married model?
#' @param add_priors_transitions Plot prior-posts for sigmas?
#' @param add_priors_trad Plot prior-posts for trad sigmas?
#' @param add_trad Plot trad parameters?
#'
#' @returns NULL, but saves a pdf with plots to the output directory
#' @keywords internal
plot_hierchecks <- function(fit, is_married, add_priors_transitions = TRUE,
                            add_priors_trad = FALSE,
                            add_trad = FALSE){
  pdf(file = file.path(fit$output_dir, "hierchecks.pdf"), width = 12, height = 8)

  if (add_priors_transitions){
    ### post and prior sigma
    for (parname in c("Omega", "Ptilde", "Betas",
                      "d_Omega", "d_Ptilde", "d_Betas")){
      print(localhierarchy::plot_prior_post_sigmas_localhierarchy(fit = fit, parname = parname))
    }
  }

  #### mu_raws
  muraw_omega1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "Omega")
  muraw_ptilde1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "Ptilde")
  muraw_splines1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "Betas", morethan1param = TRUE)
  dmuraw_omega1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "d_Omega")
  dmuraw_ptilde1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "d_Ptilde")
  dmuraw_splines1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "d_Betas", morethan1param = TRUE)

  # for Omega (1param)
  print(muraw_omega1a[["summary_plots"]][[1]])
  for (i in 1:4){
    print(muraw_omega1a[["plots_allmuraw"]][[i]])
  }
  print(dmuraw_omega1a[["summary_plots"]][[1]])
  for (i in 1:4){
    print(dmuraw_omega1a[["plots_allmuraw"]][[i]])
  }

  # Ptilde
  print(muraw_ptilde1a[["summary_plots"]][[1]])
  for (i in 1:4){
    print(muraw_ptilde1a[["plots_allmuraw"]][[i]])
  }
  print(dmuraw_ptilde1a[["summary_plots"]][[1]])
  for (i in 1:4){
    print(dmuraw_ptilde1a[["plots_allmuraw"]][[i]])
  }

  # for splines: a list with plots
  # summary plots is a list of dimension k, for each k, it gives summaries, 30 at a time
  for (k in 1:4){
    print(muraw_splines1a[["summary_plots"]][[k]][[1]])
  }
  # specific plots for one parameter are in plots[["plots_allmuraw"]]
  # muraws suggest not much updating... slightly negative in intercept
  for (k in 1:4){
    p1 = muraw_splines1a[["plots_allmuraw"]][[k]][[1]] +
      ggtitle(paste("k = ", k)) +
      theme(legend.position = "none")
    p2 = muraw_splines1a[["plots_allmuraw"]][[k]][[2]]
    p3 = muraw_splines1a[["plots_allmuraw"]][[k]][[3]]
    p4 = muraw_splines1a[["plots_allmuraw"]][[k]][[4]]
    ggall = ggarrange(p1, p2, p3, p4,
                      ncol = 2, nrow = 2,
                      common.legend = TRUE,
                      legend = "bottom")
    print(ggall)
  }
  for (k in 1:4){
    print(dmuraw_splines1a[["summary_plots"]][[k]][[1]])
  }
  # specific plots for one parameter are in plots[["plots_allmuraw"]]
  for (k in 1:4){
    p1 = dmuraw_splines1a[["plots_allmuraw"]][[k]][[1]] +
      ggtitle(paste("k = ", k)) +
      theme(legend.position = "none")
    p2 = dmuraw_splines1a[["plots_allmuraw"]][[k]][[2]]
    p3 = dmuraw_splines1a[["plots_allmuraw"]][[k]][[3]]
    p4 = dmuraw_splines1a[["plots_allmuraw"]][[k]][[4]]
    ggall = ggarrange(p1, p2, p3, p4,
                      ncol = 2, nrow = 2,
                      common.legend = TRUE,
                      legend = "bottom")
    print(ggall)
  }


  ### hierarchical parameters

  # for nond, so ds
  omega1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "Omega",
                                         hierarchical_levels = fit$hierarchical_terms_and_fixed$hierarchical_level)
  ptilde1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "Ptilde",
                                          hierarchical_levels = fit$hierarchical_terms_and_fixed$hierarchical_asymptote)
  splines1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "Betas",
                                           hierarchical_levels = fit$hierarchical_terms_and_fixed$hierarchical_splines,
                                           morethan1param =  TRUE)
  domega1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "d_Omega",
                                          hierarchical_levels = fit$d_hierarchical_terms_and_fixed$hierarchical_level)
  dptilde1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "d_Ptilde",
                                           hierarchical_levels = fit$d_hierarchical_terms_and_fixed$hierarchical_asymptote)
  dsplines1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "d_Betas",
                                            hierarchical_levels = fit$d_hierarchical_terms_and_fixed$hierarchical_splines,
                                            morethan1param =  TRUE)

  p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = omega1a)
  for (i in 1:length(p)){
    print(p[[i]] +  ggtitle(paste0("Omega, ", names(p)[i])))
  }
  p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = ptilde1a)
  for (i in 1:length(p)){
    print(p[[i]] +  ggtitle(paste0("Ptilde, ", names(p)[i])))
  }
  p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = splines1a)
  for (i in 1:length(p)){
    print(p[[i]] +  ggtitle(paste0("Betas, ", names(p)[i])))
  }
  p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = domega1a)
  for (i in 1:length(p)){
    print(p[[i]] +  ggtitle(paste0("dOmega, ", names(p)[i])))
  }
  p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = dptilde1a)
  for (i in 1:length(p)){
    p[[i]] +  ggtitle(paste0("dPtilde, ", names(p)[i]))
  }
  print(p)
  p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = dsplines1a)
  for (i in 1:length(p)){
    print(p[[i]] +  ggtitle(paste0("dBetas, ", names(p)[i])))
  }

  # for trad param
  if (add_trad){
    if (add_priors_trad ){
      ### post and prior sigma
      for (parname in c("z_Omega")){
        print(localhierarchy::plot_prior_post_sigmas_localhierarchy(fit = fit, parname = parname))
      }
    }

    #### mu_raws
    dmuraw_omega1a <- localhierarchy::plot_muraw_localhierarchy(fit = fit, parname = "z_Omega")
    print(dmuraw_omega1a[["summary_plots"]][[1]])
    for (i in 1:4){
      print(dmuraw_omega1a[["plots_allmuraw"]][[i]])
    }
    omega1a <- localhierarchy::posterior_summary_hierparam_localhierarchy(fit = fit, parname = "z_Omega",
                    hierarchical_levels = fit$z_hierarchical_terms_and_fixed$hierarchical_level)
    p <- localhierarchy::plot_posterior_summaries_localhierarchy(res = omega1a)
    for (i in 1:length(p)){
      print(p[[i]] +  ggtitle(paste0("z_Omega, ", names(p)[i])))
    }
  } # end adding trad

  dev.off()
  return(NULL)
}
