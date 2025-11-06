#' Build WorldClim predictors and a stars cube for a given country and month
#'
#' @param iso3  Character ISO3 country code (e.g., "NLD", "BEL").
#' @param month Integer month 1–12 to extract from WorldClim stacks.
#' @param vars  Character vector of variables to fetch. Default:
#'              c("tmin","tmax","prec","tavg","wind").
#' @param res   WorldClim resolution (0.5, 2.5, 5, 10). Default 0.5 arc-min.
#' @param version WorldClim version string (e.g., "2.1").
#' @param path  Download path. Default: tempdir().
#'
#' @return A list with:
#'   - stars: stars object with variables as attributes (x, y, time)
#'   - predictors: SpatRaster for the selected month, named after available vars
#'   - vars_kept: names of variables that were successfully downloaded/kept
#' @export
vsc_build_worldclim <- function(
  iso3,
  month,
  vars    = c("tmin","tmax","prec","tavg","wind"),
  res     = 0.5,
  version = "2.1",
  path    = tempdir()
) {
  stopifnot(is.character(iso3), nchar(iso3) == 3)
  if (!is.numeric(month) || length(month) != 1L || month < 1L || month > 12L) {
    stop("`month` must be a single integer from 1 to 12.")
  }

  # 1) Download each requested variable (robust to failures)
  fetch_one <- function(var) {
    tryCatch(
      geodata::worldclim_country(iso3, var = var, path = path, res = res, version = version),
      error = function(e) {
        message("[vsc_build_worldclim] Skipping variable '", var, "': ", conditionMessage(e))
        NULL
      }
    )
  }
  ras_list_raw <- lapply(vars, fetch_one)
  names(ras_list_raw) <- vars

  # 2) Keep only non-NULL with >= 12 layers (one per month)
  keep <- vapply(ras_list_raw, function(x) !is.null(x) && terra::nlyr(x) >= 12L, logical(1))
  if (!any(keep)) stop("No WorldClim variables were successfully downloaded for iso3 = ", iso3, ".")

  ras_list  <- ras_list_raw[keep]
  vars_kept <- vars[keep]

  # Ensure a clean 1..n time index for each kept stack
  for (i in seq_along(ras_list)) {
    terra::time(ras_list[[i]]) <- seq_len(terra::nlyr(ras_list[[i]]))
  }

  # 3) stars cube with ALL months (x, y, time) and variables as attributes
  stars_list  <- lapply(ras_list, stars::st_as_stars)
  stars_clima <- do.call(c, stars_list)
  stars_clima <- stats::setNames(stars_clima, vars_kept)

  # 4) SpatRaster predictors for the *selected month* (no stars slicing needed)
  #    Take the 'month'-th layer from each SpatRaster, then stack them.
  per_month <- lapply(ras_list, function(r) r[[month]])
  preds     <- terra::rast(per_month)     # layers are in the same order as vars_kept
  names(preds) <- vars_kept               # lengths now guaranteed to match

  list(
    stars      = stars_clima,   # full time cube
    predictors = preds,         # selected month
    vars_kept  = vars_kept
  )
}
