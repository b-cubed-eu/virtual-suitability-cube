#' Extract species-wise suitability for a given cell (long format)
#'
#' Works for aggregated cubes with polygon geometry (cell × species) and for
#' toy arrays without geometry (cell × species).
#'
#' @param aggregated_cube A `stars` with attribute `"suitability"` and a `species` dim.
#' @param cell_id Integer index (1-based) of the target cell.
#' @return data.frame(cell, species, suitability)
#' @export
vsc_cell_suitability_long = function(aggregated_cube, cell_id) {
  if (!inherits(aggregated_cube, "stars")) {
    stop("`aggregated_cube` must be a stars object.", call. = FALSE)
  }
  if (!("suitability" %in% names(aggregated_cube))) {
    stop("Expected attribute 'suitability'.", call. = FALSE)
  }

  dims = stars::st_dimensions(aggregated_cube)
  dn   = dimnames(aggregated_cube[["suitability"]])
  arr  = as.array(aggregated_cube[["suitability"]])

  # Trova indice specie
  if ("species" %in% names(dims)) {
    i_species = which(names(dims) == "species")
  } else {
    stop("`aggregated_cube` must have a 'species' dimension.", call. = FALSE)
  }

  # Trova indice cella (qualsiasi altra non-time/non-species) — qui assumiamo 2D cell×species
  i_other = setdiff(seq_along(dim(arr)), i_species)
  if (length(i_other) != 1L) {
    # fallback per oggetti 3D con geometria: meglio usare as.data.frame(long=TRUE)
    df = as.data.frame(aggregated_cube, long = TRUE)
    if (!all(c("species","suitability") %in% names(df))) {
      stop("Cannot reshape this stars object; expected 'species' and 'suitability' in long form.", call. = FALSE)
    }
    # crea id cella raggruppando per geometria
    if ("geometry" %in% names(df)) {
      poly_id = as.integer(factor(sf::st_as_text(df$geometry)))
      keep = which(poly_id == cell_id)
      out = df[keep, c("species","suitability")]
      out$cell = cell_id
      out = out[, c("cell","species","suitability")]
      rownames(out) = NULL
      return(out)
    }
    stop("Unsupported stars layout for this helper.", call. = FALSE)
  }

  # Estrai vettore specie per la cella richiesta
  if (i_species == 2 && i_other == 1) {
    vals = arr[cell_id, ]
  } else if (i_species == 1 && i_other == 2) {
    vals = arr[, cell_id]
  } else {
    stop("Unexpected 2D array layout; cannot locate cell/species axes.", call. = FALSE)
  }

  # Nomi specie: prova dimnames, poi st_get_dimension_values, altrimenti default
  if (!is.null(dn) && !is.null(dn[[i_species]])) {
    sp = dn[[i_species]]
  } else {
    sp = tryCatch(stars::st_get_dimension_values(aggregated_cube, "species"),
                   error = function(e) NULL)
    if (is.null(sp)) sp = paste0("sp", seq_len(length(vals)))
  }

  data.frame(cell = cell_id, species = sp, suitability = as.vector(vals), row.names = NULL)
}
