test_that("split_species_data returns named list", {
  occ <- data.frame(
    scientificName = c("sp a","sp a","sp b"),
    decimalLatitude = c(1,2,3),
    decimalLongitude = c(1,2,3)
  )
  lst <- split_species_data(occ)
  expect_type(lst, "list")
  expect_setequal(names(lst), c("sp a","sp b"))
})
