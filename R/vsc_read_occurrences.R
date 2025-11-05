#' Read GBIF occurrence file and select columns
#'
#' @param file Path to a GBIF TSV/TXT file with standard Darwin Core columns.
#' @param year_min,year_max Integer bounds to filter occurrences (inclusive).
#'
#' @return A `data.frame` with columns: scientificName, decimalLatitude, decimalLongitude, year.
#' @export
#' @importFrom utils read.delim
#' @importFrom dplyr select filter
vsc_read_occurrences = function(file, year_min = 2000, year_max = 2017) {
  occ = utils::read.delim(file, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
  keep = c("scientificName","decimalLatitude","decimalLongitude","year")
  occ = dplyr::select(occ, tidyselect::any_of(keep))
  occ = dplyr::filter(
    occ,
    !is.na(decimalLatitude), !is.na(decimalLongitude),
    !is.na(year), year >= year_min, year <= year_max
  )
  occ
}
