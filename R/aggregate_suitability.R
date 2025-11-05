#' Aggregate species suitability over a polygon grid
#'
#' Aggregates per-cell suitability for each species onto a user-provided **polygon**
#' grid (zonal statistics), computing a summary (default: mean) within each polygon.
#' This collapses the `(x, y)` spatial dimensions into a `cell` (polygon) dimension.
#'
#' @param stars_obj A `stars` object (e.g., from `vsc_predict_sdm_for_new_area`) with
#'   a suitability attribute and a `species` dimension. Attribute must be `"suit"`
#'   or `"suitability"`.
#' @param grid An `sf` polygon layer (POLYGON/MULTIPOLYGON) defining zones.
#' @param fun Summary function applied within polygons (default `mean`).
#' @param na_rm Logical; remove NAs inside `fun` (default `TRUE`).
#' @return A `stars` object with attribute `"suitability"`, polygon-based spatial
#'   dimension (from `grid`) and a `species` dimension.
#' @export
#' @importFrom stars st_dimensions st_redimension st_set_dimensions st_as_stars st_get_dimension_values
#' @importFrom sf st_geometry_type st_crs st_transform st_bbox st_as_sfc
#' @importFrom terra rast extract values<- crs<- nlyr ext
aggregate_suitability <- function(stars_obj, grid, fun = mean, na_rm = TRUE) {
  if (!inherits(stars_obj, "stars")) {
    stop("`stars_obj` must be a stars object.", call. = FALSE)
  }
  if (!inherits(grid, "sf")) {
    stop("`grid` must be an sf polygon layer (POLYGON/MULTIPOLYGON).", call. = FALSE)
  }
  gtype <- unique(as.character(sf::st_geometry_type(grid, by_geometry = TRUE)))
  if (!all(gtype %in% c("POLYGON","MULTIPOLYGON"))) {
    stop("`grid` must be POLYGON/MULTIPOLYGON.", call. = FALSE)
  }

  # pick suitability attribute
  attr_name <- intersect(c("suit","suitability"), names(stars_obj))
  if (!length(attr_name)) stop("Expected suitability attribute 'suit' or 'suitability'.", call. = FALSE)
  attr_name <- attr_name[1]

  # species dimension
  if (!("species" %in% names(stars::st_dimensions(stars_obj)))) {
    stop("`stars_obj` must have a 'species' dimension.", call. = FALSE)
  }
  species_vals <- tryCatch(
    stars::st_get_dimension_values(stars_obj, "species"),
    error = function(e) NULL
  )
  if (is.null(species_vals)) {
    d <- stars::st_dimensions(stars_obj)
    nsp <- d[["species"]]$to - d[["species"]]$from + 1
    species_vals <- paste0("sp", seq_len(nsp))
  }

  # armonizza CRS
  crs_obj <- sf::st_crs(stars_obj)
  if (!is.na(crs_obj)) {
    if (is.na(sf::st_crs(grid))) {
      sf::st_crs(grid) <- crs_obj
    } else if (sf::st_crs(grid) != crs_obj) {
      grid <- sf::st_transform(grid, crs_obj)
    }
  }

  # tieni solo l'attributo di suitability e splitta per specie
  suit_only  <- stars_obj[attr_name]
  by_species <- split(suit_only, "species")  # list of stars, one per species

  # helper: converti stars (1 attr) -> SpatRaster, ricostruendo se serve
  .stars_to_rast <- function(x) {
    rr <- try(terra::rast(x), silent = TRUE)
    if (!inherits(rr, "try-error")) return(rr)

    # ricostruzione manuale da bbox + dims
    d  <- stars::st_dimensions(x)
    bb <- sf::st_bbox(x)
    # individua le prime due dimensioni spaziali
    dim_names <- names(d)
    if (length(d) < 2) stop("stars object lacks 2 spatial dimensions.", call. = FALSE)
    nx <- d[[1]]$to - d[[1]]$from + 1
    ny <- d[[2]]$to - d[[2]]$from + 1
    rr <- terra::rast(ncol = nx, nrow = ny,
                      xmin = bb["xmin"], xmax = bb["xmax"],
                      ymin = bb["ymin"], ymax = bb["ymax"])
    # valori: prendi il primo attributo
    vals <- as.vector(x[[1]])
    terra::values(rr) <- vals
    if (!is.na(sf::st_crs(x))) terra::crs(rr) <- sf::st_crs(x)
    rr
  }

  # zonal stats via terra::extract per ogni specie
  agg_list <- vector("list", length(by_species))
  for (i in seq_along(by_species)) {
    sp_st <- by_species[[i]]
    rr    <- .stars_to_rast(sp_st)
    vals  <- terra::extract(rr, grid, fun = fun, na.rm = na_rm)
    vals  <- vals[order(vals$ID), , drop = FALSE]
    sf_i  <- sf::st_sf(suitability = vals[[2]], geometry = sf::st_geometry(grid))
    agg_i <- stars::st_as_stars(sf_i)   # stars con geometria poligonale
    agg_list[[i]] <- agg_i
  }
  names(agg_list) <- as.character(species_vals)

  combined <- do.call(c, agg_list)      # più attributi (uno per specie)
  combined <- stars::st_redimension(combined)
  combined <- stars::st_set_dimensions(
    combined, which = "new_dim", values = names(agg_list), names = "species"
  )
  names(combined) <- "suitability"
  combined
}
