#' Download WorldClim and build predictors (stars + SpatRaster for a month)
#'
#' Downloads monthly WorldClim variables for a country code (ISO3), converts them
#' to a `stars` object with attributes `tmin`, `tmax`, `prec`, `tavg`, `wind`,
#' and extracts a training month as `SpatRaster` predictors.
#'
#' @param iso3 Character ISO3 country code (e.g., "NLD" or "Nld").
#' @param variables Character vector of WorldClim variables. Defaults to
#'   `c("tmin","tmax","prec","tavg","wind")`. Valid values include
#'   `tmin`, `tmax`, `tavg`, `prec`, `wind`, `vapr`, `bio`.
#' @param res Numeric WorldClim resolution in arc-minutes. One of 10, 5, 2.5, 0.5.
#' @param version Character WorldClim version, e.g. "2.1".
#' @param month Integer in 1:12; which month to extract as training predictors.
#' @param path Character output directory for `geodata::worldclim_country()` (temporary by default).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{stars_climate}: `stars` object with 3 dims (x,y,time) and attributes = `variables`.
#'   \item \code{predictors}: `SpatRaster` with the selected `month` and layer names = `variables`.
#' }
#'
#' @examples
#' \dontrun{
#' cl <- vsc_build_worldclim("NLD", month = 5)
#' cl$predictors
#' }
#' @export
#' @importFrom geodata worldclim_country
#' @importFrom stars st_as_stars
#' @importFrom terra rast
vsc_build_worldclim <- function(
  iso3,
  variables = c("tmin","tmax","prec","tavg","wind"),
  res = 0.5,
  version = "2.1",
  month = 5,
  path = tempdir()
) {
  stopifnot(is.character(iso3), length(iso3) == 1L)
  stopifnot(month %in% 1:12)

  # download each variable as SpatRaster (12 layers each)
  rasters <- lapply(variables, function(v) {
    geodata::worldclim_country(iso3, v, path = path, res = res, version = version)
  })

  # convert each to stars, bind as attributes and name them
  stars_list <- lapply(rasters, stars::st_as_stars)
  stars_clima <- do.call(c, stars_list)
  names(stars_clima) <- variables

  # month
  clima_m <- stars_clima[, , month, drop = FALSE]
  predictors <- terra::rast(clima_m)
  names(predictors) <- variables

  list(stars_climate = stars_clima, predictors = predictors)
}
