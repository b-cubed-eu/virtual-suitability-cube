test_that("aggregate_suitability aggregates over polygon grid", {
  r1 <- terra::rast(ncol = 4, nrow = 4, xmin = 3, xmax = 5, ymin = 50, ymax = 52, vals = runif(16))
  r2 <- terra::rast(ncol = 4, nrow = 4, xmin = 3, xmax = 5, ymin = 50, ymax = 52, vals = runif(16))
  terra::crs(r1) <- "EPSG:4326"
  terra::crs(r2) <- "EPSG:4326"

  s1 <- stars::st_as_stars(r1)
  s2 <- stars::st_as_stars(r2)
  st <- c(s1, s2)                         # due attributi
  st <- stars::st_redimension(st)         # nuova dimensione "new_dim"
  st <- stars::st_set_dimensions(st, which = "new_dim",
                                 values = c("sp1", "sp2"), names = "species")
  names(st) <- "suit"                     # attributo di suitability

  bb  <- sf::st_as_sfc(sf::st_bbox(st))
  pg  <- sf::st_make_grid(bb, n = c(2, 2), what = "polygons")
  g   <- sf::st_sf(geometry = pg)
  sf::st_crs(g) <- sf::st_crs(st)
  g$id <- seq_len(nrow(g))

  out <- vscube::aggregate_suitability(st, g, fun = mean)

  expect_s3_class(out, "stars")
  expect_equal(names(out), "suitability")
  expect_true("species" %in% names(stars::st_dimensions(out)))
})
