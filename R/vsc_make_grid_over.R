#' Build a polygon grid over the extent of an object
#'
#' Creates an `sf` polygon grid over the bounding box of `x`. You can define
#' the grid by a target cellsize or by number of cells per side.
#'
#' @param x A `SpatRaster`, `stars`, or `sf` object providing the extent/CRS.
#' @param cellsize Numeric length (map units) of grid cells. If `NULL`, use `n`.
#' @param n Integer vector of length 2: number of cells in x and y (used when
#'   `cellsize` is `NULL`). Default `c(50, 50)`.
#' @param square Logical; if `FALSE`, allow non-square cells if needed.
#'
#' @return An `sf` polygon grid with an `id` column.
#'
#' @examples
#' \dontrun{
#' cl_bel <- vsc_build_worldclim("BEL", month = 5)
#' grd <- vsc_make_grid_over(cl_bel$predictors, n = c(50, 50))
#' }
#' @export
#' @importFrom sf st_bbox st_as_sfc st_crs st_make_grid st_sf
#' @importFrom stars st_as_stars
vsc_make_grid_over <- function(x, cellsize = NULL, n = c(50, 50), square = TRUE) {
  # bbox + crs robusti per diversi tipi
  bb  <- sf::st_bbox(x)
  sfc <- sf::st_as_sfc(bb)
  crs <- sf::st_crs(x)

  pg <- if (!is.null(cellsize)) {
    sf::st_make_grid(sfc, cellsize = cellsize, what = "polygons", square = square)
  } else {
    sf::st_make_grid(sfc, n = n,                what = "polygons", square = square)
  }
  g <- sf::st_sf(geometry = pg)
  sf::st_crs(g) <- crs
  g$id <- seq_len(nrow(g))
  g
}
