test_that("vsc_cell_id_for_point returns an id inside grid", {
  r  <- terra::rast(ncol=2, nrow=2, vals=1:4)
  bb <- sf::st_as_sfc(sf::st_bbox(r))
  pg <- sf::st_make_grid(bb, n = c(2,2))
  g  <- sf::st_sf(geometry = pg); sf::st_crs(g) <- 4326
  g$id <- seq_len(nrow(g))

  id <- vsc_cell_id_for_point(g, lon = 0, lat = 0)  # dentro bbox (se bbox è long/lat standard)
  expect_true(is.integer(id) || is.numeric(id))
  expect_true(!is.na(id))
})

test_that("vsc_cell_suitability_long returns long df for one cell", {
  arr <- array(1:4, dim = c(cell=1, species=4))
  s   <- stars::st_as_stars(arr); names(s) <- "suitability"
  s   <- stars::st_set_dimensions(s, "species", values = paste0("sp", 1:4))

  df <- vsc_cell_suitability_long(s, 1L)
  expect_true(all(c("cell","species","suitability") %in% names(df)))
  expect_equal(nrow(df), 4)
  expect_equal(df$cell[1], 1)
})
