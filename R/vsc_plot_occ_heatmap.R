#' Plot heatmap of occurrences per species and year
#'
#' @param occ A data.frame from \code{vsc_read_occurrences()}.
#' @param break_by Integer step for legend breaks (default auto).
#' @return A `ggplot` object.
#' @export
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn labs theme_minimal theme element_text coord_fixed
#' @importFrom viridisLite mako
vsc_plot_occ_heatmap <- function(occ, break_by = NULL) {
  stopifnot(all(c("scientificName","year") %in% names(occ)))
  counts <- dplyr::ungroup(
    dplyr::summarise(
      dplyr::group_by(occ, scientificName, year),
      count = dplyr::n(),
      .groups = "drop"
    )
  )
  maxc <- if (nrow(counts)) max(counts$count, na.rm = TRUE) else 0
  breaks <- if (is.null(break_by) && maxc > 0) seq(0, maxc, length.out = min(6, maxc + 1)) else seq(0, maxc, by = break_by %||% 1)
  ggplot2::ggplot(counts, ggplot2::aes(x = year, y = scientificName, fill = count)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = viridisLite::mako(10, direction = -1), breaks = breaks, labels = breaks) +
    ggplot2::labs(title = "Occurrences per species per year", x = "year", y = "species", fill = "occurrences") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "right",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::coord_fixed(ratio = 1.5)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
