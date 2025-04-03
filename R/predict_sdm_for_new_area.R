#' Predicting suitability in a new area for the same species using the models trained previously 
#'
#' This function creates, for a geographic area defined by the user, a new SDM map for the same species
#' using new predictions and the previous models 
#' 
#' @param models A list of MaxEnt models, one for each species. T
#' hey must be the output of create_sdm_for_species_list function
#' @param stack_clima A SpatRaster object with environmental variables for the new area

#' @return A stars R object with suitability as attribute and x, y, and species as dimension


predict_sdm_for_new_area <- function(models, new_stack) {
  predictions <- list()
  
  for (species_name in names(models)) {
    message("SDM for: ", species_name)
    
    # Estrai il modello per la specie
    mx <- models[[species_name]]
    
    # Applica il modello al nuovo stack
    pred_map <- predictEnmSdm(mx, new_stack)
    
    # Rinomina il layer come "suitability"
    names(pred_map) <- "suitability"
    
    # Salva la previsione nella lista
    predictions[[species_name]] <- pred_map
  }
  
  # Converte la lista di SpatRaster in un unico SpatRaster
  prediction_stack <- rast(predictions)
  
  # Converte in oggetto stars con dimensioni corrette
  stars_predictions <- st_as_stars(prediction_stack) %>% 
    setNames("suit") %>% 
    st_set_dimensions(., names = c("x", "y", "species"))
  
  return(stars_predictions)
}
