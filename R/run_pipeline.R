#' Run the vscube SDM pipeline
#'
#' Download climate, read & split occurrences, train MaxEnt for all species, and predict suitability.
#'
#' @param occ_file Path to GBIF TSV/TXT (scientificName, decimalLatitude, decimalLongitude, year).
#' @param iso3 ISO3 country code for climate (e.g., "NLD").
#' @param month Integer 1–12: which month to use as predictors (default 5 = May).
#' @param background_points Integer background size for MaxEnt (default 10000).
#' @param variables Character vector of WorldClim variables (default: c("tmin","tmax","prec","tavg","wind")).
#' @param res,version,path Passed to \code{vsc_build_worldclim()}.
#' @param predictors Optional character vector: predictor names (default all).
#' @param verbose Logical; print progress (default TRUE).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{climate}: list(stars_climate, predictors)
#'   \item \code{occ}: cleaned occurrence table
#'   \item \code{species_list}: split occurrences by species (list of data.frames)
#'   \item \code{sdms}: list(models, predictions)
#' }
#' @export
vsc_run_pipeline = function(
  occ_file,
  iso3,
  month = 5,
  background_points = 10000,
  variables = c("tmin","tmax","prec","tavg","wind"),
  res = 0.5,
  version = "2.1",
  path = tempdir(),
  predictors = NULL,
  verbose = TRUE
) {
  cl  = vsc_build_worldclim(iso3, variables = variables, res = res, version = version, month = month, path = path)
  occ = vsc_read_occurrences(occ_file)
  sp  = split_species_data(occ)
  sdms = vsc_create_sdm_for_species_list(
    sp, cl$predictors,
    background_points = background_points,
    predictors = predictors,
    verbose = verbose
  )
  list(climate = cl, occ = occ, species_list = sp, sdms = sdms)
}
