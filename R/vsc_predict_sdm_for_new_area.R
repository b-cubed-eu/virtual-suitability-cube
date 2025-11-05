#' Predict SDMs for a new area using previously trained models
#'
#' Applies a list of fitted models (one per species) to a new `SpatRaster`
#' of predictors, returning a `stars` cube with a `species` dimension
#' and a single `suit` attribute.
#'
#' @param models Named list of fitted models (e.g., from
#'   \code{vsc_create_sdm_for_species_list()}), one per species.
#' @param new_stack A `SpatRaster` with environmental predictors for the new area.
#' @param predict_fun Function used to perform the prediction. Default:
#'   \code{enmSdmX::predictEnmSdm}. It must accept arguments (model, new_stack)
#'   and return a `SpatRaster`.
#' @param verbose Logical; print progress messages. Default `TRUE`.
#'
#' @return A `stars` object with attribute `"suit"` and dimensions `x`, `y`, `species`.
#' @export
#' @importFrom stars st_as_stars st_set_dimensions
#' @importFrom terra rast nlyr ncol nrow ext
vsc_predict_sdm_for_new_area <- function(models,
                                         new_stack,
                                         predict_fun = enmSdmX::predictEnmSdm,
                                         verbose = TRUE) {
  if (!inherits(new_stack, "SpatRaster")) {
    stop("`new_stack` must be a SpatRaster.", call. = FALSE)
  }
  if (!length(models)) {
    stop("`models` is empty.", call. = FALSE)
  }
  if (!is.function(predict_fun)) {
    stop("`predict_fun` must be a function(model, new_stack) -> SpatRaster.", call. = FALSE)
  }

  preds <- vector("list", length(models))
  nms   <- names(models)
  if (is.null(nms) || any(!nzchar(nms))) nms <- paste0("sp", seq_along(models))

  for (i in seq_along(models)) {
    sp <- nms[i]
    if (isTRUE(verbose)) message("[vscube] Predicting: ", sp)
    mx <- models[[i]]

    r  <- predict_fun(mx, new_stack)  # << injection point
    if (!inherits(r, "SpatRaster")) {
      stop("`predict_fun` must return a SpatRaster.", call. = FALSE)
    }
    names(r) <- "suitability"
    preds[[i]] <- r
  }
  names(preds) <- nms

  stk <- terra::rast(preds)              # stack dei layer (uno per specie)
  s   <- stars::st_as_stars(stk)
  names(s) <- "suit"
  s   <- stars::st_set_dimensions(s, names = c("x", "y", "species"))
  s
}
