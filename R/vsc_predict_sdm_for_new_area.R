#' Project trained maxnet models to a new area
#'
#' Applies a list of trained **maxnet** models (one per species) to a new
#' predictor stack (SpatRaster) and returns a `stars` cube with attribute
#' `"suit"` and dimension `"species"`.
#'
#' @param models    Named list of fitted maxnet models, e.g. the output of
#'   [vsc_create_sdm_for_species_list()]$models. Each element name should be
#'   the species name.
#' @param new_stack A `terra::SpatRaster` with the predictor layers required by
#'   the models. Layer names should match the variables used in training.
#'
#' @return A `stars` object with attribute `"suit"` and dimensions `x`, `y`,
#'   and `species`. One layer per species.
#' @export
vsc_predict_sdm_for_new_area <- function(models, new_stack) {
  # ---- sanity checks ---------------------------------------------------------
  if (length(models) == 0L) {
    stop("No models provided (length(models) == 0). Train at least one model first.")
  }
  if (!inherits(new_stack, "SpatRaster")) {
    stop("`new_stack` must be a terra::SpatRaster with predictor layers.")
  }

  # ---- predict per-species ---------------------------------------------------
  pred_list <- list()

  for (sp in names(models)) {
    message("[vscube] Predicting: ", sp)
    mx <- models[[sp]]

    pr <- try(
      terra::predict(
        object = new_stack,
        model  = mx,
        fun    = function(m, v) {
          # v arrives as matrix/data.frame of cell covariates
          # maxnet exposes predict.maxnet via stats::predict()
          as.numeric(stats::predict(m, as.data.frame(v), type = "cloglog"))
        },
        na.rm  = TRUE
      ),
      silent = TRUE
    )

    if (inherits(pr, "try-error")) {
      warning("[vscube] Skipping '", sp, "' during prediction: ",
              conditionMessage(attr(pr, "condition")))
      next
    }

    names(pr) <- "suitability"
    pred_list[[sp]] <- pr
  }

  if (length(pred_list) == 0L) {
    stop("All species failed during projection (no predictions to combine). ",
         "Check that new_stack has the required predictor names.")
  }

  # ---- combine to stars cube -------------------------------------------------
  # Combine SpatRasters (one per species) into multilayer SpatRaster
  rr <- terra::rast(pred_list)           # layers named by species
  st <- stars::st_as_stars(rr)           # -> stars
  st <- stats::setNames(st, "suit")      # attribute name

  # label dims as x, y, species (third dim is often "band")
  dims <- names(stars::st_dimensions(st))
  if (length(dims) == 3L) {
    st <- stars::st_set_dimensions(st, names = c("x", "y", "species"))
  }
  st
}
