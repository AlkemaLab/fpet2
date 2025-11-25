#' Standardize country names
#'
#' @param name_country country name
#'
#' @returns Standardized country name
#' @keywords internal
std_name_country <- function(name_country) {
  name_country <- gsub(" ", "_", name_country, fixed = TRUE)
  name_country <- gsub("(", "", name_country, fixed = TRUE)
  name_country <- gsub(")", "", name_country, fixed = TRUE)
  name_country <- gsub(",", "", name_country, fixed = TRUE)
  name_country <- gsub("'", "", name_country, fixed = TRUE)
  name_country <- gsub("-", "", name_country, fixed = TRUE)
  name_country <- gsub("ç", "c", name_country, fixed = TRUE)
  name_country <- gsub("È", "E", name_country, fixed = TRUE)
  name_country <- gsub("ô", "o", name_country, fixed = TRUE)

  return(name_country)
}
