#' Quick raster + grid plot (ggplot)
#'
#' @param r A `SpatRaster` (single-layer) to plot as background.
#' @param grid An `sf` polygon grid.
#' @param layer Optional layer name or index if `r` has multiple layers.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_raster geom_sf scale_fill_viridis_c theme_void theme element_blank margin
#' @importFrom terra as.data.frame
vsc_plot_raster_with_grid = function(r, grid, layer = 1) {
  if (!inherits(r, "SpatRaster")) stop("`r` must be a SpatRaster.", call. = FALSE)
  if (nlyr(r) > 1) r = r[[layer]]

  df = terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  val_col = names(df)[3L]  # nome colonna valori

  ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = .data[[val_col]])) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::geom_sf(data = grid, color = "black", size = 0.4, fill = NA, alpha = 0.2) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
