#' Train MaxEnt models and predict suitability for a list of species
#'
#' For each species in `species_list`, trains a MaxEnt model (`enmSdmX::trainMaxEnt`)
#' using presence-background data derived from `stack_clima`, then predicts suitability.
#'
#' @param species_list A named list of data.frames (one per species) with columns
#'   `decimalLatitude` and `decimalLongitude` in lon/lat (EPSG:4326).
#' @param stack_clima A `SpatRaster` with environmental predictors (same extent/CRS used for all species).
#' @param background_points Integer number of background samples per species (default 10000).
#' @param predictors Optional character vector: layer names to use as predictors. Default: all layers.
#' @param verbose Logical; if `TRUE`, prints progress messages.
#'
#' @return A list with two named lists:
#'   \itemize{
#'     \item \code{models}: fitted MaxEnt models per species.
#'     \item \code{predictions}: `SpatRaster` predictions per species.
#'   }
#'
#' @examples
#' \dontrun{
#' cl = vsc_build_worldclim("NLD", month = 5)
#' occ = vsc_read_occurrences("occurrence.txt")
#' sp  = split_species_data(occ)
#' res = vsc_create_sdm_for_species_list(sp, cl$predictors, background_points = 5000)
#' }
#' @export
#' @importFrom terra vect crds cellFromXY extract spatSample
#' @importFrom enmSdmX trainMaxEnt predictEnmSdm elimCellDuplicates
vsc_create_sdm_for_species_list = function(
  species_list,
  stack_clima,
  background_points = 10000,
  predictors = NULL,
  verbose = TRUE
) {
  if (!inherits(stack_clima, "SpatRaster")) {
    stop("`stack_clima` must be a SpatRaster.", call. = FALSE)
  }
  if (is.null(predictors)) predictors = names(stack_clima)

  models = list()
  predictions = list()

  for (species_name in names(species_list)) {
    if (isTRUE(verbose)) message("[vscube] Processing species: ", species_name)

    dat = species_list[[species_name]]
    # drop rows without coords
    dat = dat[!is.na(dat$decimalLatitude) & !is.na(dat$decimalLongitude), , drop = FALSE]
    if (!nrow(dat)) {
      if (isTRUE(verbose)) message("[vscube] Skipping ", species_name, " (no valid coordinates).")
      next
    }

    # create points as SpatVector in lon/lat
    pts = terra::vect(dat, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

    # keep only points that overlay predictors cells
    idx = !is.na(terra::cellFromXY(stack_clima, terra::crds(pts)))
    pts = pts[idx, ]
    if (nrow(pts) < 5L) {  # arbitrary small threshold
      if (isTRUE(verbose)) message("[vscube] Skipping ", species_name, " (too few points after filtering).")
      next
    }

    # remove duplicated cells (one presence per cell)
    pts = enmSdmX::elimCellDuplicates(pts, stack_clima)

    # extract env at presences
    occEnv = terra::extract(stack_clima, pts, ID = FALSE)
    occEnv = occEnv[stats::complete.cases(occEnv), , drop = FALSE]

    # sample background and extract env
    bg = terra::spatSample(stack_clima, background_points)
    bg   = bg[stats::complete.cases(bg), , drop = FALSE]
    if (nrow(bg) > background_points) bg = bg[seq_len(background_points), , drop = FALSE]

    # response + bind env
    presBg = data.frame(presBg = c(rep(1, nrow(occEnv)), rep(0, nrow(bg))))
    env    = rbind(occEnv, bg)
    env    = cbind(presBg, env)

    # ensure requested predictors exist
    missing_pred = setdiff(predictors, names(env))
    if (length(missing_pred)) {
      stop("Missing predictors in data: ", paste(missing_pred, collapse = ", "), call. = FALSE)
    }

    # train MaxEnt
    mx = enmSdmX::trainMaxEnt(
      data = env,
      resp = "presBg",
      predictors = predictors,
      verbose = isTRUE(verbose)
    )
    models[[species_name]] = mx

    # predict
    mxMap = enmSdmX::predictEnmSdm(mx, stack_clima)
    predictions[[species_name]] = mxMap

    if (isTRUE(verbose)) message("[vscube] Done: ", species_name)
  }

  list(models = models, predictions = predictions)
}
