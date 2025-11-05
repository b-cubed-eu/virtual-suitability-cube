test_that("vsc_predict_sdm_for_new_area works with injected predict_fun", {
  # SpatRaster minimale (3x3, 1 layer)
  r <- terra::rast(ncol = 3, nrow = 3, vals = runif(9))
  names(r) <- "bio1"

  fake_models <- list(a = list(), b = list())

  fake_predict <- function(model, new_stack) {
    terra::rast(ncol = terra::ncol(new_stack),
                nrow = terra::nrow(new_stack),
                ext  = terra::ext(new_stack),
                vals = runif(terra::ncol(new_stack) * terra::nrow(new_stack)))
  }

  out <- vsc_predict_sdm_for_new_area(fake_models, r, predict_fun = fake_predict, verbose = FALSE)

  expect_s3_class(out, "stars")
  expect_equal(names(out), "suit")
  expect_true("species" %in% names(stars::st_dimensions(out)))
})
