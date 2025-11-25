



#' Write stan models
#'
#' @param add_aggregates Add national aggregates? makes sense only for subnational all-region run for 1 country
#' @param add_trad Include trad use estimation and data?
#' @param all_women Estimate for all_women or by marital group?
#' @param add_emu Include EMU estimation and data?
#'
#' @returns writes a stan file to the stan directory. The model files are named
#' based on the options chosen.
#' @export
#'

# testing: for no aggr, add_trad
write_model <- function(add_aggregates = FALSE,
                        add_trad = TRUE,
                        all_women = FALSE,
                        add_emu = FALSE){

  # general
  # make sure that names are not used in function calls etc!
  # and part of other things... eg eta = bad!! ie in beta
  # we use capitalization to avoid issues

  ### build base model
  stan_code <- readr::read_file(here::here("inst/stan", "fpem_buildingblocks.stan"))

  ### add information on the observations
  data_data_code <- readr::read_file(here::here("inst/stan", "data_data.stan"))
  data_data_allwomen_code <- readr::read_file(here::here("inst/stan", "data_data_allwomen.stan"))
  stan_code <- stan_code %>%
      stringr::str_replace_all("\\{\\{DATA_DATA\\}\\}",
                               ifelse(all_women,
                                      data_data_allwomen_code,
                                      data_data_code))

  # for all remaining blocks, they need to be duplicated for unmarried, with unmarried_ added to relevant parameters


  ### building blocks for process models demand and demand satisfied
  processmodel <- list(
    processmodel_data_code = readr::read_file(here::here("inst/stan", "processmodel_data.stan")),
    processmodel_transformeddata_code = readr::read_file(here::here("inst/stan", "processmodel_transformeddata.stan")),
    processmodel_parameters_code = readr::read_file(here::here("inst/stan", "processmodel_parameters.stan")),
    processmodel_transformedparameters_code = readr::read_file(here::here("inst/stan", "processmodel_transformedparameters.stan")),
    processmodel_model_code = readr::read_file(here::here("inst/stan", "processmodel_model.stan")),
    processmodel_generatedquantities_code = readr::read_file(here::here("inst/stan", "processmodel_generatedquantities.stan"))
  )
  processmodel_demand <- processmodel
  processmodel_tmp <- processmodel
  params <- c("Epsilon", "Eta", "Betas", "Ptilde", "Omega", "Rho", "Tau")
  for (block in names(processmodel)){
    processmodel_demand[[block]] <- replacesomething_inblock(processmodel[[block]], params, "d_")
    processmodel_tmp[[block]] <- replacesomething_inblock(processmodel[[block]], params, "tmp_")
  }
  if (all_women){
    processmodel_demand <- map(processmodel_demand, function(x) add_unmarried_inblock(x))
    processmodel_tmp <- map(processmodel_tmp, function(x) add_unmarried_inblock(x))
  }
  # plug it in
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DEMAND_DATA\\}\\}",
                             processmodel_demand$processmodel_data_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DEMAND_TRANSFORMEDDATA\\}\\}",
                             processmodel_demand$processmodel_transformeddata_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DEMAND_PARAMETERS\\}\\}",
                             processmodel_demand$processmodel_parameters_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DEMAND_TRANSFORMEDPARAMETERS\\}\\}",
                             processmodel_demand$processmodel_transformedparameters_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DEMAND_MODEL\\}\\}",
                             processmodel_demand$processmodel_model_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DEMAND_GENERATEDQUANTITIES\\}\\}",
                             processmodel_demand$processmodel_generatedquantities_code)
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_DATA\\}\\}",
                             processmodel_tmp$processmodel_data_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRANSFORMEDDATA\\}\\}",
                             processmodel_tmp$processmodel_transformeddata_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_PARAMETERS\\}\\}",
                             processmodel_tmp$processmodel_parameters_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRANSFORMEDPARAMETERS\\}\\}",
                             processmodel_tmp$processmodel_transformedparameters_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_MODEL\\}\\}",
                             processmodel_tmp$processmodel_model_code) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_GENERATEDQUANTITIES\\}\\}",
                             processmodel_tmp$processmodel_generatedquantities_code)

  #### add in data and data models parameters
  datamodel_data_code <- readr::read_file(here::here("inst/stan", "datamodel_data.stan"))
  datamodel_parameters_code <- readr::read_file(here::here("inst/stan", "datamodel_parameters.stan"))
  datamodel_transformedparameters_code <- readr::read_file(here::here("inst/stan", "datamodel_transformedparameters.stan"))
  datamodel_model_param_code <- readr::read_file(here::here("inst/stan", "datamodel_model_param.stan"))
  #datamodel_generatedquantities_code <- readr::read_file(here::here("inst/stan", "datamodel_generatedquantities.stan"))
  if (all_women){
    datamodel_data_code <- add_unmarried_inblock(datamodel_data_code)
    datamodel_parameters_code <- add_unmarried_inblock(datamodel_parameters_code)
    datamodel_transformedparameters_code <- add_unmarried_inblock(datamodel_transformedparameters_code)
    datamodel_model_param_code <- add_unmarried_inblock(datamodel_model_param_code)
  }
  # start with the 2s
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{DATAMODEL2_DATA\\}\\}", datamodel_data_code) %>%
    stringr::str_replace_all("\\{\\{DATAMODEL2_PARAMETERS\\}\\}", datamodel_parameters_code) %>%
    stringr::str_replace_all("\\{\\{DATAMODEL2_TRANSFORMEDPARAMETERS\\}\\}", datamodel_transformedparameters_code)%>%
    stringr::str_replace_all("\\{\\{DATAMODEL2_MODEL_PARAM\\}\\}", datamodel_model_param_code)
  # replace DM1 by DM2
  stan_code <- stan_code %>%
    stringr::str_replace_all("DM1", "DM2")

  # now consider dm3, depends on whether it is used
  if (!add_trad){
    stan_code <- stan_code %>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_DATA\\}\\}", " ") %>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_PARAMETERS\\}\\}", " ") %>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_TRANSFORMEDPARAMETERS\\}\\}", " ")%>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_MODEL_PARAM\\}\\}", " ")
  } else {
    stan_code <- stan_code %>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_DATA\\}\\}", datamodel_data_code) %>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_PARAMETERS\\}\\}", datamodel_parameters_code) %>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_TRANSFORMEDPARAMETERS\\}\\}", datamodel_transformedparameters_code)%>%
      stringr::str_replace_all("\\{\\{DATAMODEL3_MODEL_PARAM\\}\\}", datamodel_model_param_code)
    # replace DM1 by DM3
    stan_code <- stan_code %>%
      stringr::str_replace_all("DM1", "DM3")
  }

  # now add blocks for DM1
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{DATAMODEL_DATA\\}\\}", datamodel_data_code) %>%
    stringr::str_replace_all("\\{\\{DATAMODEL_PARAMETERS\\}\\}", datamodel_parameters_code) %>%
    stringr::str_replace_all("\\{\\{DATAMODEL_TRANSFORMEDPARAMETERS\\}\\}", datamodel_transformedparameters_code) %>%
    stringr::str_replace_all("\\{\\{DATAMODEL_MODEL_PARAM\\}\\}", datamodel_model_param_code)

  # the data model itself depends on whether trad is added and if it's for all_women and if there are aggregates
  datamodel_model_code <- readr::read_file(here::here("inst/stan", "datamodel_model.stan"))
  datamodel_model_wtrad_code <- readr::read_file(here::here("inst/stan", "datamodel_model_wtrad.stan"))
  datamodel_model_wtrad_code_allwomen <- readr::read_file(here::here("inst/stan", "datamodel_model_wtrad_allwomen.stan"))
  datamodel_model_wtrad_code_allwomen_noaggregates <-
    readr::read_file(here::here("inst/stan", "datamodel_model_wtrad_allwomen_noaggregates.stan"))

  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{DATAMODEL_MODEL\\}\\}",
             ifelse(all_women,
                    # for all women
                    ifelse(add_aggregates, datamodel_model_wtrad_code_allwomen, datamodel_model_wtrad_code_allwomen_noaggregates),
                    # not all women
                    ifelse(add_trad, datamodel_model_wtrad_code, datamodel_model_code)))


  ### consider aggregates
  aggregates_data_code <- readr::read_file(here::here("inst/stan", "aggregates_data.stan"))
  aggregates_generatedquantities_code <- readr::read_file(here::here("inst/stan", "aggregates_generatedquantities.stan"))
  if (all_women){
    aggregates_data_code <- add_unmarried_inblock(aggregates_data_code)
    aggregates_generatedquantities_code <- add_unmarried_inblock(aggregates_generatedquantities_code)
  }
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{AGGREGATES_DATA\\}\\}",
                             ifelse(add_aggregates, aggregates_data_code, " ")) %>%
    stringr::str_replace_all("\\{\\{AGGREGATES_GENERATEDQUANTITIES\\}\\}",
                             ifelse(add_aggregates, aggregates_generatedquantities_code, " "))

  ### trad processmodel related info
  # already did the data and data model parameters and data model above
  # now the process model
  processmodel_trad_data_code <- readr::read_file(here::here("inst/stan", "processmodel_trad_data.stan"))
  processmodel_trad_transformeddata_code <- readr::read_file(here::here("inst/stan", "processmodel_trad_transformeddata.stan"))
  processmodel_trad_parameters_code <- readr::read_file(here::here("inst/stan", "processmodel_trad_parameters.stan"))
  processmodel_trad_transformedparameters_code <- readr::read_file(here::here("inst/stan", "processmodel_trad_transformedparameters.stan"))
  processmodel_trad_model_code <- readr::read_file(here::here("inst/stan", "processmodel_trad_model.stan"))
  if (all_women){
    processmodel_trad_data_code <- add_unmarried_inblock(processmodel_trad_data_code)
    processmodel_trad_transformeddata_code <- add_unmarried_inblock(processmodel_trad_transformeddata_code)
    processmodel_trad_parameters_code <- add_unmarried_inblock(processmodel_trad_parameters_code)
    processmodel_trad_transformedparameters_code <- add_unmarried_inblock(processmodel_trad_transformedparameters_code)
    processmodel_trad_model_code <- add_unmarried_inblock(processmodel_trad_model_code)
  }
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRAD_DATA\\}\\}",
                             ifelse(add_trad, processmodel_trad_data_code, " ")) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRAD_TRANSFORMEDDATA\\}\\}",
                             ifelse(add_trad, processmodel_trad_transformeddata_code, " ")) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRAD_PARAMETERS\\}\\}",
                             ifelse(add_trad, processmodel_trad_parameters_code, " ")) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRAD_TRANSFORMEDPARAMETERS\\}\\}",
                             ifelse(add_trad, processmodel_trad_transformedparameters_code, " ")) %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRAD_MODEL\\}\\}",
                             ifelse(add_trad, processmodel_trad_model_code, " "))

  # generated quantiies trad depends on if trad is added and if it's for all_women
  processmodel_trad_generatedquantities_code <- readr::read_file(here::here("inst/stan", "processmodel_trad_generatedquantities.stan"))
  processmodel_trad_generatedquantities_code_allwomen <-
                readr::read_file(here::here("inst/stan", "processmodel_trad_generatedquantities_allwomen.stan"))
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{PROCESSMODEL_TRAD_GENERATEDQUANTITIES\\}\\}",
                             ifelse(all_women,
                                    processmodel_trad_generatedquantities_code_allwomen,
                                    ifelse(add_trad, processmodel_trad_generatedquantities_code,
                                           " ")))

  ### all_women code blocks
  allwomen_generatedquantities_code_allwomen <-
    readr::read_file(here::here("inst/stan", "allwomen_generatedquantities.stan"))
  aggregatesandallwomen_generatedquantities_code_allwomen <-
    readr::read_file(here::here("inst/stan", "aggregatesandallwomen_generatedquantities.stan"))

  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{ALLWOMEN_GENERATEDQUANTITIES\\}\\}",
                             ifelse(all_women,
                                    allwomen_generatedquantities_code_allwomen, " ")) %>%
    stringr::str_replace_all("\\{\\{AGGREGATESALLWOMEN_GENERATEDQUANTITIES\\}\\}",
                             ifelse(all_women & add_aggregates,
                                    aggregatesandallwomen_generatedquantities_code_allwomen, " "))


  ## remove "tmp_"
  stan_code <- stan_code %>%
    stringr::str_replace_all("tmp_", "")

  # emus?
  emu_data_code <- readr::read_file(here::here("inst/stan", "emu_data.stan"))
  emu_parameters_code <- readr::read_file(here::here("inst/stan", "emu_parameters.stan"))
  emu_model_code <- readr::read_file(here::here("inst/stan", "emu_model.stan"))
  stan_code <- stan_code %>%
    stringr::str_replace_all("\\{\\{EMU_DATA\\}\\}",
                             ifelse(add_emu, emu_data_code, " ")) %>%
    stringr::str_replace_all("\\{\\{EMU_PARAMETERS\\}\\}",
                             ifelse(add_emu, emu_parameters_code, " ")) %>%
    stringr::str_replace_all("\\{\\{EMU_MODEL\\}\\}",
                             ifelse(add_emu, emu_model_code, " "))

  # write to file
  # all_women true only if add_trad is true
  if (! add_trad & !add_aggregates){
    complete_model <- cmdstanr::write_stan_file(stan_code, dir = here::here("inst/stan"),
                                                force_overwrite = TRUE,
                                               # basename = "fpem_buildingblocks_complete")
                                               basename = "fpem_notallwomen_notrad")
  }
  if (add_trad & !all_women){
    complete_model <- cmdstanr::write_stan_file(stan_code, dir = here::here("inst/stan"),
                                                force_overwrite = TRUE,
                                                # basename = "fpem_buildingblocks_complete")
                                                basename = "fpem_notallwomen_wtrad")
  }
  if (add_trad & all_women & !add_aggregates & !add_emu){
    complete_model <- cmdstanr::write_stan_file(stan_code, dir = here::here("inst/stan"),
                                                force_overwrite = TRUE,
                                                #basename = "fpem_buildingblocks_complete_wtrad_allwomen")
                                                basename = "fpem_allwomen")
  }
  if (add_trad & all_women & !add_aggregates & add_emu){
    complete_model <- cmdstanr::write_stan_file(stan_code, dir = here::here("inst/stan"),
                                                force_overwrite = TRUE,
                                                #basename = "fpem_buildingblocks_complete_wtrad_allwomen_emu")
                                                basename = "fpem_allwomen_emu")
  }
  if (all_women & add_aggregates  & !add_emu){
    complete_model <- cmdstanr::write_stan_file(stan_code, dir = here::here("inst/stan"),
                                                force_overwrite = TRUE,
                                                #basename = "fpem_buildingblocks_complete_allwomenaggregates")
                                                basename = "fpem_allwomen_aggregates")

  }
  if (all_women & add_aggregates & add_emu){
    complete_model <- cmdstanr::write_stan_file(stan_code, dir = here::here("inst/stan"),
                                                force_overwrite = TRUE,
                                                #basename = "fpem_buildingblocks_complete_allwomenaggregates_emu")
                                                basename = "fpem_allwomen_aggregates_emu")
  }
  return(invisible(NULL))
}

#' Replace something in block
#'
#' @param block string with stan code block
#' @param params vector of parameter names to replace
#' @param prefix prefix to add to parameter names
#'
#' @returns modified block
#'
#' @keywords internal
replacesomething_inblock <- function(block, params = c("Epsilon", "Eta",
                                               "Betas", "Ptilde", "Omega", "Rho", "Tau"),
                             prefix = "d_"){
  # make sure that names are not used in function calls etc!
  # and part of other things... eg eta = bad!! ie in beta
  # we use capitalization to avoid issues
  prefixed_param_name <- paste0(prefix, params)
  for (i in 1:length(prefixed_param_name)){
    block <- block %>%
      stringr::str_replace_all(params[i], prefixed_param_name[i])
  }
  return(block)
}

#' Get all parameters to replace
#'
#' @returns vector of parameter names to replace
#' @keywords internal
get_all_params_toreplace <- function(){
  params_to_collect <- dplyr::bind_rows(
    expand.grid(
      prefix = c("tmp_", "d_", "z_"),
      param = c("Omega", "Eta", "Epsilon",
                "Rho", "Tau")),
    expand.grid(
      prefix = c("tmp_", "d_", "z_"),
      param = c("Betas", "Ptilde"))
    #,
    #expand.grid(
    #    prefix = c("DM1_", "DM2_", "DM3_"),
    #    param = c("nonse", "sdbias", "rho_pma", "local_shrinkage_dm", "global_shrinkage_dm", "caux_dm"))
  ) |>
    dplyr::mutate(prefixed_param_name = paste0(prefix, param)) %>%
    pull(prefixed_param_name)
  return(params_to_collect)
}




add_unmarried_inblock <- function(block,
                                  unmarried_prefix = "unmarried_"){
  paramstoreplace <- c(get_all_params_toreplace(), "z_fix_smoothing",
                       "N", "S", "isdhs", "source",
                       c("DM1_", "DM2_", "DM3_"),
                       c("geo_unit_natpop_weight_tr",
                       #  "demand_aggr", "modern_aggr", "unmet_aggr", "trad_aggr"))
                       "demand", "modern", "unmet", "trad"))
  block_unmarried <- block
  for (i in 1:length(paramstoreplace)){
    block_unmarried <- block_unmarried %>%
      stringr::str_replace_all(paramstoreplace[i],
                               paste0(unmarried_prefix, paramstoreplace[i]))
  }
  paste(block, block_unmarried, sep = "\n")
}
#cmdstanr::write_stan_file(paste("bla", "bla2", sep = "\n"), here())

# write_model(add_aggregates = FALSE,
#             add_trad = TRUE,
#             all_women = TRUE)
# #fpem_buildingblocks_complete_wtrad_allwomen
# write_model(add_aggregates = FALSE,
#             add_trad = TRUE,
#             all_women = FALSE)
