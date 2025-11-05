test_that("vsc_read_occurrences selects and filters columns", {
  tf <- tempfile(fileext = ".tsv")
  df <- data.frame(
    scientificName = c("A", "B"),
    decimalLatitude = c(10, NA),
    decimalLongitude = c(20, 30),
    year = c(2005, 1995),
    extra = 1:2
  )
  write.table(df, tf, sep = "\t", row.names = FALSE, quote = FALSE)
  out <- vsc_read_occurrences(tf, year_min = 2000, year_max = 2010)
  expect_true(all(c("scientificName","decimalLatitude","decimalLongitude","year") %in% names(out)))
  expect_true(all(out$year >= 2000 & out$year <= 2010))
  expect_false(any(is.na(out$decimalLatitude)))
})
