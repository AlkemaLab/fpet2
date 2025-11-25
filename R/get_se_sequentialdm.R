
#' Get effective sample size from proportion and its standard error
#'
#' @param p proportion
#' @param se_p standard error of the proportion
#'
#' @returns effective sample size
#' @keywords internal
get_effective_samplesize <- function(p, se_p){
  1/se_p^2*p*(1-p)
}

#' Get Binomial standard error from proportion and sample size
#'
#' @param p proportion
#' @param n sample size
#'
#' @returns standard error
#' @keywords internal
get_se <- function(p, n){
  sqrt(p*(1-p)/n)
}


#' Get standard error on logit scale from proportion and its standard error
#'
#' @param prop proportion
#' @param se_prop  standard error of the proportion (original scale)
#'
#' @returns standard error on logit scale
#' @keywords internal
get_se_logitprop <- function(prop, se_prop){
  se_prop/(prop - prop^2)
}

#' probit <- function(x) pnorm(x)
#' # stan
#' #phiapprox <- function(x) inv_logit(0.07056 * x ^ 3 + 1.5976 * x)
#' inv_probit <- function(x) qnorm(x)
#'
#' #' Get standard error on inverse probit scale from proportion and its standard error
#' #'
#' #' @param prop proportion
#' #' @param se_prop standard error of the proportion (original scale)
#' #'
#' #' @returns standard error on inverse probit scale
#' #' @keywords internal
#' get_se_invprobitprop <- function(prop, se_prop){
#'   # se_prop*abs(derivate evaluated at prop)
#'   se_prop*1/dnorm(inv_probit(prop))
#' }

# logit_wbounds <- function(x, lower, upper) {
#   log((x-lower)/(upper - x))
# }
#
# inv_logit_wbounds <- function(x, lower, upper) {
#   lower + (upper - lower)/(1+exp(-x))
# }

#p <- seq(0.05, 0.95, 0.01)
#invlogit_wbounds(logit_wbounds(p, 0.01, 0.99), 0.01, 0.99)

#' Logit function
#'
#' @param x  input value
#'
#' @returns logit of x
#' @keywords internal
logit <- function(x) {
  log(x/(1-x))
}

#' Inverse logit function
#'
#' @param x input value
#'
#' @returns inverse logit of x
#' @keywords internal
inv_logit <- function(x) {
  1/(1+exp(-x))
}



#' Get standard errors for FP indicators in sequential DM approach
#'
#' @param nat_data data frame with survey data
#' @param use_small_SE_setting logical, if TRUE (default), use smaller effective sample sizes when imputing SEs
#'
#' @returns tibble with SEs added
#' @keywords internal
get_se_sequentialdm <- function(nat_data,
                                use_small_SE_setting = TRUE
                                ){
  nat_data <- nat_data %>%
    rename(se_unmet_any = se_unmet_need) %>%
    mutate(
      se_modern = ifelse(se_modern == 0, NA, se_modern),
      se_unmet_any = ifelse(se_unmet_any == 0, NA, se_unmet_any),
      se_traditional = ifelse(se_traditional == 0, NA, se_traditional))

  # get effective samples sizes, to later use to impute missing SEs
  nat_data <- nat_data %>%
    mutate(neff_mod = get_effective_samplesize(contraceptive_use_modern, se_modern),
           neff_trad = get_effective_samplesize(contraceptive_use_traditional, se_traditional),
           neff_unmet_need_any = get_effective_samplesize(unmet_any, se_unmet_any))

  if ( use_small_SE_setting){
    nat_data <- nat_data %>%
      mutate( neff_derived = pmax(neff_mod, neff_trad, neff_unmet_need_any, na.rm = TRUE))
  } else {
    nat_data <- nat_data %>%
      mutate( neff_derived = pmin(neff_mod, neff_trad, neff_unmet_need_any, na.rm = TRUE))
  }

  # approach:
  # for dhs, use effective sample size from most recent survey with nonNA sample size
  print("When imputing SEs for DHS, we use the effective sample size of its preceding survey")
  # use record_id+fixed to merge back in later
  neff_dhs <-  nat_data %>%
    dplyr::filter(data_series_type == "DHS") %>%
    dplyr::mutate(neff_dhs = neff_derived) %>%
    dplyr::group_by(name_country, region_code) %>%
    dplyr::mutate(no = n()) %>%
    dplyr::filter(no > 1) %>%
    dplyr::arrange(start_date) %>%
    dplyr::mutate(
      neff_dhs_prior = dplyr::lag(neff_dhs, 1),  # still NA for first DHS but that's fine
      # also NA if more than 1 is missing but not going to worry about that now either
      neff_dhs_use = ifelse(is.na(neff_dhs), neff_dhs_prior, neff_dhs)
    ) %>%
    ungroup() %>%
    dplyr::select(record_id_fixed, neff_dhs_prior, neff_dhs, neff_dhs_use)

  # neff_country <-  nat_data %>%
  #   group_by(name_country, region_code) %>%
  #   mutate(neff_country = min(neff_derived)) %>%
  #   ungroup() %>%
  #   dplyr::select(record_id_fixed, neff_country)

  neff_country <-  nat_data %>%
    dplyr::mutate(neff_country = ifelse(region_code %in% c("National", "") | is.na(region_code), 500, 200)) %>%
    dplyr::select(record_id_fixed, neff_country)

  # no longer needed, kept just in case
  # for populations without any info, impute a fixed number for effective sample size
  if ( use_small_SE_setting){
   # print("When imputing SEs for populations without any SEs, we assume an effective sample size of 500")
    neff_impute = 500
  } else {
   # print("When imputing SEs for countries without any SEs, we assume an effective sample size of 200")
    neff_impute = 200
  }

  nat_data <- nat_data %>%
    dplyr::left_join(neff_dhs) %>%
    dplyr::left_join(neff_country) %>%
    dplyr::mutate(
      # effective sample sizes
      neff_use = case_when(
        # first check if we have observed something for other indicators
        !is.na(neff_derived) ~ neff_derived,
        # then DHS
        data_series_type == "DHS" & !is.na(neff_dhs_use) ~ neff_dhs_use,
        # then country level
          !is.na(neff_country) ~ neff_country,
        # then impute
        TRUE ~ neff_impute),
      neff_denominator_1minmcpr = neff_use*(1-contraceptive_use_modern),
      neff_denominator_1minuse = neff_use*(1-use_any),
      # ses
      se_modern = ifelse(!is.na(se_modern), se_modern,
                         get_se(p = contraceptive_use_modern, n = neff_use)),
      se_unmet_any = ifelse(!is.na(se_unmet_any), se_unmet_any,
                         get_se(p = unmet_any, n = neff_use)),
      se_traditional = ifelse(!is.na(se_traditional), se_traditional,
                            get_se(p = contraceptive_use_traditional, n = neff_use)),
      # for the ratios
      se_unmet_overnotmodern =  get_se(p = unmet_overnotmodern,
                                        n = neff_denominator_1minmcpr),
      se_trad_overnotmodern =  get_se(p = trad_overnotmodern,
                                       n = neff_denominator_1minmcpr),
      se_trad_overunmet =  get_se(p = trad_overunmet,
                                         n = neff_use*unmet),
      se_unmetany_overnonuse =  get_se(p = unmetany_overnonuse,
                                      n = neff_denominator_1minuse),
       # get logits
      logit_contraceptive_use_modern = logit(contraceptive_use_modern),
      se_logit_modern = get_se_logitprop(contraceptive_use_modern, se_modern),
      logit_unmet_overnotmodern = logit(unmet_overnotmodern),
      se_logit_unmet_overnotmodern =  get_se_logitprop(unmet_overnotmodern, se_unmet_overnotmodern),
      logit_trad_overnotmodern = logit(trad_overnotmodern),
      se_logit_trad_overnotmodern =  get_se_logitprop(trad_overnotmodern, se_trad_overnotmodern),
      logit_trad_overunmet = logit(trad_overunmet),
      se_logit_trad_overunmet =  get_se_logitprop(trad_overunmet,se_trad_overunmet),
      logit_unmetany_overnonuse = logit(unmetany_overnonuse),
      se_logit_unmetany_overnonuse = get_se_logitprop(unmetany_overnonuse, se_unmetany_overnonuse)
       )

  # set se's to a minimum value
  se_minimum <- 0.01
  apply_minsd <- function(sds)
    ifelse(sds < se_minimum, se_minimum, sds)
  nat_data <- nat_data %>%
    dplyr::mutate(
           se_logit_modern = apply_minsd(se_logit_modern),
           se_logit_unmet_overnotmodern =  apply_minsd( se_logit_unmet_overnotmodern),
           se_logit_trad_overnotmodern = apply_minsd(se_logit_trad_overnotmodern),
           se_logit_trad_overunmet =  apply_minsd(se_logit_trad_overunmet),
           se_logit_unmetany_overnonuse =  apply_minsd(se_logit_unmetany_overnonuse)
    )
  nat_data
}

