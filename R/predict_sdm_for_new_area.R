#' Predicting suitability in a new area for the same species using the models trained previously 
#'
#' This function creates, for a geographic area defined by the user, a new SDM map for the same species
#' using new predictions and the previous models 
#' 
#' @param models A list of MaxEnt models, one for each species. T
#' hey must be the output of `create_sdm_for_species_list` function
#' @param stack_clima A SpatRaster object with environmental variables for the new area

#' @return A stars R object with suitability as attribute and x, y, and species as dimension


predict_sdm_for_new_area = function(models, new_stack) {
  predictions = list()
  
  for (species_name in names(models)) {
    message("SDM for: ", species_name)
    
    # extract model for species 
    mx = models[[species_name]]
    
    # apply to new stack
    pred_map = predictEnmSdm(mx, new_stack)
    
    # rename layer as 'suitability'
    names(pred_map) = "suitability"
    
    # save
    predictions[[species_name]] = pred_map
  }
  
  # from list to stack
  prediction_stack = rast(predictions)
  
  # from stack to rast
  stars_predictions = st_as_stars(prediction_stack) %>% 
    setNames("suit") %>% 
    st_set_dimensions(., names = c("x", "y", "species"))
  
  return(stars_predictions)
}
