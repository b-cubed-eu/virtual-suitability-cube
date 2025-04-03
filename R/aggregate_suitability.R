#' Aggregating suitability values of the species over polygons for a given area
#'
#' This function aggregates in a grid defined by the user the suitability values of the species
#' contained in the data cube previously generated with `predict_sdm_for_new_area` function. 
#' Within each polygon, an average of the suitability of each species is calculated. This way, the (x,y) dimensions 
#' are collapsed into a cell dimension. 
#' 
#' @param stars_obj A stars object generated with `predict_sdm_for_new_area` that contains suitability values for each species
#' in the new area 
#' @param grid A polygon over which aggregate suitability values by mean

#' @return A stars R object with suitability as attribute and cells and species as dimension

aggregate_suitability <- function(stars_obj, grid) {
  
  # extract species name
  species_names <- st_get_dimension_values(stars_obj, "species")
  
  agg_list <- vector("list", length(species_names))
  
  for (i in seq_along(species_names)) {
    species_name <- species_names[i]
    message("Aggregating suitability for: ", species_name)
    
    # extract suitability for the species
    species_raster <- stars_obj["suit", , , i, drop = TRUE]
    
    # aggregate over grid by mean
    agg_suit <- aggregate(species_raster, grid, mean, as_points = TRUE, na.action = na.omit)
    
    # save
    agg_list[[i]] <- agg_suit
  }
  
  # combine everything into a stars object
  aggregated_cube <- do.call(c, agg_list) %>%
    st_redimension() %>%
    st_set_dimensions(., which = "new_dim", values = species_names, names = "species") %>%
    setNames("suitability")
  
  return(aggregated_cube)
}
