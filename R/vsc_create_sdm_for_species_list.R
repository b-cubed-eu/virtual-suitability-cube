#' Train MaxEnt-like SDMs for each species (pure R maxnet backend)
#'
#' @param species_list list of data.frames, one per species (from split_species_data)
#' @param stack_clima  SpatRaster of predictors for the training area/month
#' @param background_points integer number of background points
#' @param predictors character vector of predictor layer names (defaults to all)
#' @param verbose logical
#' @return list(models = <list of maxnet models>, predictions = <list of SpatRaster>)
#' @export
vsc_create_sdm_for_species_list <- function(
  species_list,
  stack_clima,
  background_points = 10000,
  predictors = NULL,
  verbose = TRUE
) {
  if (is.null(predictors)) predictors <- names(stack_clima)

  models <- list()
  predictions <- list()

  for (species_name in names(species_list)) {
    if (verbose) message("[vscube] Processing species: ", species_name)
    dat <- species_list[[species_name]]

    # keep only rows with coords
    dat <- dat[!is.na(dat$decimalLatitude) & !is.na(dat$decimalLongitude), ]
    if (nrow(dat) < 2L) {
      if (verbose) message("  - too few occurrences, skipping")
      next
    }

    # presence points as SpatVector (WGS84)
    pres <- terra::vect(dat,
                        geom = c("decimalLongitude", "decimalLatitude"),
                        crs  = "EPSG:4326")

    # drop points falling outside predictor cells
    pres <- pres[!is.na(terra::cellFromXY(stack_clima, terra::crds(pres))), ]

    # drop duplicate presences within same raster cell  << FIXED HERE
    cells <- terra::cellFromXY(stack_clima, terra::crds(pres))
    keep  <- !duplicated(cells)
    pres  <- pres[keep, ]

    # extract env at presences
    occEnv <- terra::extract(stack_clima, pres, ID = FALSE)
    occEnv <- occEnv[complete.cases(occEnv), , drop = FALSE]
    if (nrow(occEnv) < 2L) {
      if (verbose) message("  - too few valid env rows, skipping")
      next
    }

    # background
    bgEnv <- terra::spatSample(stack_clima, size = background_points, method = "random")
    bgEnv <- bgEnv[complete.cases(bgEnv), , drop = FALSE]
    if (nrow(bgEnv) > background_points) bgEnv <- bgEnv[seq_len(background_points), , drop = FALSE]

    # assemble pres/bg
    presBg <- c(rep(1L, nrow(occEnv)), rep(0L, nrow(bgEnv)))
    envDF  <- rbind(occEnv, bgEnv)
    envDF  <- envDF[, predictors, drop = FALSE]

    # fit maxnet (pure R)
    # default feature classes often ok; you can expose 'f' or 'regmult' if needed
    fm <- maxnet::maxnet.formula(p = presBg, data = envDF, classes = "lqph")
    mx <- maxnet::maxnet(p = presBg, data = envDF, f = fm)  # no Java at all

    models[[species_name]] <- mx

    # in-area prediction for QA (training geography)
    pred_map <- terra::predict(
      object = stack_clima,
      model  = mx,
      fun    = function(m, v) {
        # v is a matrix/data.frame of cell values -> coerce to data.frame
        as.numeric(stats::predict(m, as.data.frame(v), type = "cloglog"))
      },
      na.rm  = TRUE
    )


    names(pred_map) <- "suitability"
    predictions[[species_name]] <- pred_map
  }

  list(models = models, predictions = predictions)
}
