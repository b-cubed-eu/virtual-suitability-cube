test_that("vsc_create_sdm_for_species_list runs on toy data", {
  testthat::skip_if_not_installed("enmSdmX")
  # toy 3x3 predictors
  r <- terra::rast(ncol=3, nrow=3, vals=runif(9))
  r <- c(r, terra::rast(ncol=3, nrow=3, vals=runif(9))) # 2 layers
  names(r) <- c("p1","p2")

  sp <- list(
    "sp1" = data.frame(decimalLatitude = c(0.1, 0.2, 0.3), decimalLongitude = c(0.1, 0.2, 0.3))
  )
  res <- vsc_create_sdm_for_species_list(sp, r, background_points = 50, predictors = names(r), verbose = FALSE)
  expect_true(is.list(res$models))
  expect_true(is.list(res$predictions))
})
