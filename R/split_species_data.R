#' Split a dataset by species
#'
#' This function takes a dataframe containing multiple species and
#' divides it by species
#'
#' @param occ_data An object of class `data.frame`, it must contain  
#' scientificName column, the standard format given by gbif
#' @return A list of n dataframes, one per species

split_species_data = function(occ_data) {
  
  # check if the dataset has 'scientificName' column
  if(!"scientificName" %in% colnames(occ_data)) {
    stop("Dataset does not have scientificName column'.")
  }
  
  # species list
  species_list = split(occ_data, occ_data$scientificName)
  
  # iteration
  for (species in names(species_list)) {
    species_name = gsub(" ", "_", species)  # replace spaces with "_"
    assign(species_name, species_list[[species]], envir = .GlobalEnv)
  }
  
  message("Done")
  
  # species list
  return(species_list)
}
