#' Plot per-cell suitability per species (scatter)
#'
#' @param df A data.frame with columns `species` and `suitability` (e.g., from
#'   `vsc_cell_suitability_long()`).
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme element_text
vsc_plot_cell_suitability = function(df) {
  if (!all(c("species","suitability") %in% names(df))) {
    stop("`df` must have columns 'species' and 'suitability'.", call. = FALSE)
  }
  ggplot2::ggplot(df, ggplot2::aes(x = species, y = suitability, color = species)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12)
    )
}
