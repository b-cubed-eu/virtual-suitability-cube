test_that("vsc_build_worldclim structure looks right", {
  testthat::skip("Skipping download in CI")
  cl <- vsc_build_worldclim("NLD", month = 5)
  expect_s3_class(cl$stars_climate, "stars")
  expect_s4_class(cl$predictors, "SpatRaster")
})
