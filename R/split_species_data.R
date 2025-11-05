#' Split a GBIF-like occurrence table into a list by species
#'
#' @param occ_data A `data.frame` with at least `scientificName`, `decimalLatitude`, `decimalLongitude`.
#' @return A named `list(data.frame, ...)` with one element per species.
#' @export
split_species_data =function(occ_data) {
  if (!"scientificName" %in% names(occ_data)) {
    stop("`occ_data` must contain a 'scientificName' column.", call. = FALSE)
  }
  split(occ_data, occ_data$scientificName)
}
