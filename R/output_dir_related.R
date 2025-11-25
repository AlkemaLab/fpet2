
#' Get directory of output, located in parent of wd
#'
#' @param ... character string subdirectories of the bayestransition_output
#' folder
#'
#' @returns character string of the full path to the output directory
#' @keywords internal
#'
get_output_dir <- function(...) {
  do.call(file.path,
          c(list(here::here() %>% dirname(), "bayestransition_output"), ...))
}



#' Get relative output directory path
#'
#' Constructs a relative path to a subfolder inside `bayestransition_output`.
#'
#' @param folder_name Character. The name of the subfolder inside `bayestransition_output`.
#'
#' @returns character string of the full path to the output directory
#' @keywords internal
get_relative_output_dir <- function(folder_name) {
  file.path("..", "bayestransition_output", folder_name)
}

#' Rename output folder and directory in fit object
#'
#' Renames an existing folder in `bayestransition_output`and updates the `output_dir` in the model fit object.
#'
#' @param indicator Character. Indicator name.
#' @param run_step Character. Run step - either 1a, 1b, local_national.
#' @param old_folder_name Character. The current folder name (inside `bayestransition_output`).
#' @param new_folder_name Character. The new name to rename the folder to.
#'
#' @returns NULL. The function performs the renaming and updates the fit object in place.
#'
#' @keywords internal
rename_output_folder <- function(indicator, run_step, old_folder_name, new_folder_name) {

  old_dir <- get_relative_output_dir(old_folder_name)
  new_dir <- get_relative_output_dir(new_folder_name)

  if (!dir.exists(old_dir)) stop("Old directory doesn't exist: ", old_dir)
  if (dir.exists(new_dir)) stop("New directory already exists: ", new_dir)

  success <- file.rename(old_dir, new_dir)
  if (!success) stop("Failed to rename the folder.")

  # load and update fit object
  fit_path <- file.path(new_dir, paste0(indicator, "_fit_wpostsumm.rds"))
  fit <- readRDS(fit_path)
  fit$output_dir <- new_dir
  saveRDS(fit, fit_path)

  # if summary exists, update
  summary_path <- file.path(new_dir, paste0(indicator, "_summary", run_step, ".rds"))
  if (file.exists(summary_path)) {
    summary <- readRDS(summary_path)
    summary$output_dir <- new_dir
    saveRDS(summary, summary_path)
  }
}

#' Load fit
#'
#' Load a fit object. This fit should be saved in a folder that contains
#' a file named `fit.rds` with the fpemplus object as well as csv files produced
#' by cmdstanr as a result of the estimation process.
#'
#' @param path path to the folder containing the model fit
#'
#' @return object of class fpemplus
#'
#' @keywords internal
load_fit <- function(path) {
  fit <- readRDS(Sys.glob(file.path(path, "fit*.rds"))[1]) # just take first fit object
  sample_files <- Sys.glob(file.path(path, "*.csv"))
  sample_files <- sample_files[!grepl("diagnostic", sample_files)]
  fit$samples <- cmdstanr::as_cmdstan_fit(sample_files)

  return(fit)
}


#' Load a saved model fit
#'
#' This helper function loads a previously saved model fit based on the indicator name,
#' run step, and folder suffix.
#'
#' @param indicator Character string specifying the indicator name (e.g., `"ideliv"`, `"anc4"`).
#' @param runstep Character string indicating the model step (e.g., `"step1a"`, `"step1b"`).
#' @param folder_suffix Character or numeric. The suffix at the end of the folder name.
#'
#' @return A model fit object read from an RDS file.
#' @export
get_fit <- function(indicator, runstep, folder_suffix) {
  runname <- paste0(indicator, "_", runstep, "_", folder_suffix)
  file_path <- file.path(get_output_dir(runname), paste0(indicator, "_fit_wpostsumm.rds"))

  if (!file.exists(file_path)) {
    stop("The specified fit file does not exist: ", file_path)
  }

  fit <- readRDS(file_path)
  fit$output_dir <- file.path(here::here() %>% dirname(), "bayestransition_output", runname)

  return(fit)
}

