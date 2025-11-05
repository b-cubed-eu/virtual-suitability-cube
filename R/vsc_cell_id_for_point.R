#' Get grid cell ID for a lon/lat point
#'
#' Returns the `id` (row index or an `id` column if present) of the polygon
#' in `grid` that contains the point (`lon`, `lat` in EPSG:4326).
#'
#' @param grid An `sf` polygon layer (POLYGON/MULTIPOLYGON).
#' @param lon,lat Numeric longitude/latitude (EPSG:4326).
#' @param id_col Optional name of the id column in `grid`. If `NULL`, returns the row index.
#'
#' @return A single integer id (or `NA_integer_` if no polygon contains the point).
#' @export
#' @importFrom sf st_sfc st_point st_as_sfc st_set_crs st_crs st_transform st_join
vsc_cell_id_for_point = function(grid, lon, lat, id_col = "id") {
  if (!inherits(grid, "sf")) stop("`grid` must be an sf layer.", call. = FALSE)
  pt = sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  # proietta al CRS della griglia se necessario
  if (!is.na(sf::st_crs(grid)) && sf::st_crs(grid)$epsg != 4326) {
    pt = sf::st_transform(pt, sf::st_crs(grid))
  }
  hit = suppressWarnings(sf::st_join(sf::st_as_sf(pt), grid))
  if (!nrow(hit)) return(NA_integer_)
  if (!is.null(id_col) && id_col %in% names(hit)) {
    as.integer(hit[[id_col]][1])
  } else {
    as.integer(attr(hit, "row.names")[1])
  }
}
