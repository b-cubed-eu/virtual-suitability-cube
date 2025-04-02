#' Creating SDMs with MaxEnt for many species in the same area
#'
#' This function creates, for a geographic area defined by the user, a MaxEnt model
#' that can be implemented to make predictions. 
#' The function takes as input a list of species (output from the previous function)
#' a stack of predictors, i. e. climatic layers for the same area. MaxEnt comes from
#' `enmSdmX` package. 
#' 
#' @param species A list of species. They must have decimalLatitude and decimalLongitude as columns
#' @param stack_clima A SpatRaster object with environmental variables
#' @param background_points The number of background points that the user wants to generate within MaxEnt
#' @param predictors A string of names of the predictors that will be used. If not specified, all the environemntal
#' variables in the stack will be used as predictors
#' coordinates as c(xmin, xmax, ymin, ymax)
#' @return an object with $models and $predictions for each species 


create_sdm_for_species_list <- function(species, stack_clima, background_points = 10000, predictors = NULL) {
  
  # empty list
  models <- list()
  predictions <- list()
  
  # loop for each species
  for (species_name in names(species)) {
    message("Processing species: ", species_name)
    
    species_data <- species[[species_name]]
    
    # NA remove
    species_data <- species_data[!is.na(species_data$decimalLatitude) & !is.na(species_data$decimalLongitude), ]
    
    # SpatVector object with coordinates
    species_spat <- vect(species_data, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    
    # clean points outside cells
    species_spat <- species_spat[!is.na(terra::cellFromXY(stack_clima, terra::crds(species_spat))), ]
    
    # delete duplicates
    species_spat <- elimCellDuplicates(species_spat, stack_clima)
    
    # values extraction
    occEnv <- terra::extract(stack_clima, species_spat, ID = FALSE)
    occEnv <- occEnv[complete.cases(occEnv), ]
    
    # background points
    bgEnv <- terra::spatSample(stack_clima, background_points)
    bgEnv <- bgEnv[complete.cases(bgEnv), ]
    bgEnv <- bgEnv[1:min(background_points, nrow(bgEnv)), ]
    
    # presence-background in the same df
    presBg <- data.frame(
      presBg = c(
        rep(1, nrow(occEnv)),
        rep(0, nrow(bgEnv))
      )
    )
    
    # combine
    env <- rbind(occEnv, bgEnv)
    env <- cbind(presBg, env)
    
    if (is.null(predictors)) {
      predictors <- names(stack_clima)
    }
    
    #  MaxEnt model from enmSdmX
    mx <- trainMaxEnt(
      data = env,
      resp = 'presBg',
      predictors = predictors,
      verbose = TRUE
    )
    
    # save in list 'models'
    models[[species_name]] <- mx
    
    # prediction for the species with enmSdmX
    mxMap <- predictEnmSdm(mx, stack_clima)
    
    # save map prediction in 'predictions'
    predictions[[species_name]] <- mxMap
    
    message("done for: ", species_name)  # Messaggio intermedio al termine della specie
  }
  
  return(list(models = models, predictions = predictions))
}
