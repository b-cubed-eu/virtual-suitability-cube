# load packages
library(sf)
library(raster)
library(viridis)
library(virtualspecies)
library(ggplot2)
library(tidyverse)
library(terra)
library(stars)
library(geodata)
library(purrr)
library(reshape2)

## Download climatic data 
# worldclim_country function. 
# The res paramater is the resolution, with valid values as 10,5,2.5 and 0.5 
# The var parameter is the variable name, with valid values as tmin, tmax, tavg, prec, wind, vapr and bio.
# The worldclim_country function will return 12 layers of raster (in a SpatRaster), for each variable. 

tmin <- worldclim_country("Lux", "tmin", path=tempdir(), res = 0.5, version = "2.1")
tmax <- worldclim_country("Lux", "tmax", path=tempdir(), res = 0.5, version = "2.1")
prec <- worldclim_country("Lux", "prec", path=tempdir(), res = 0.5, version = "2.1")

# 12 months of precipitation in Luxembourg. Each layer represents the average value for a month.
plot(tmin, col = plasma(500, alpha = 1, begin = 0, end = 1, direction = 1))

# here, we will simply use numbers from 1 to 12 to represent # the months.
time(tmin) <- time(tmax) <- time(prec) <- 1:12
tmin

climate_vars <- c("tmin", "tmax", "prec")

# print(tmin)

## Stars Data Cube
# first, we create a list that contains all the stacks, to which we apply st_as_stars.
stars_clima <- list(tmin, tmax, prec) %>%
  lapply(st_as_stars) %>%
  do.call("c", .) %>%
  setNames(climate_vars)

print(stars_clima)

## Virtual Species
# Let's choose to work with the first two months. For both months, we will create a stack containing the 3 climatic variables.

# number of months
num_months <- 2


# list of RasterStacks for each month
raster_stacks <- map(1:num_months, function(month) {
  
  # extraction and convertion of the climatic variables
  layers <- map(climate_vars, function(var) {
    as(stars_clima[var,,,month], "Raster")  
    # convertion into rasters
  })
  
  # layer stack of the current month
  stack_month <- stack(layers)
  # names
  names(stack_month) <- climate_vars
  return(stack_month)
})


january <- raster_stacks[[1]]
february <- raster_stacks[[2]]

print(january)

generate_suitability <- function(climate_stack, seed_value) {
  set.seed(seed_value)
  random_sp <- generateRandomSp(raster.stack = climate_stack,
                                convert.to.PA = FALSE,
                                species.type = "multiplicative",
                                approach = "response",
                                relations = "gaussian",
                                realistic.sp = TRUE,
                                plot = FALSE)
  
  return(random_sp$suitab.raster)
}

# generate the suitability for species 1 in January and February and create a single object for species 1 that contains the suitability for both months
suit_sp1_jan <- generate_suitability(january, seed_value = 121)
suit_sp1_feb <- generate_suitability(february, seed_value = 121)
suit_sp1 <- c(suit_sp1_feb, suit_sp1_jan)  

suit_sp1

# same for species 2
suit_sp2_jan <- generate_suitability(january, seed_value = 456)
suit_sp2_feb <- generate_suitability(february, seed_value = 456)
suit_sp2 <- c(suit_sp2_feb, suit_sp2_jan)  

# add time dimension
time(suit_sp2) <- time(suit_sp1) <- 1:2  

# suitability for Species 1
plot(suit_sp1, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
title("Suitability of Species 1 in Jan-Feb", outer=TRUE, line=-0.9)

# suitability for Species 2
plot(suit_sp2, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
title("Suitability of Species 2 in Jan-Feb", outer=TRUE, line=-0.9)

# stars cube for Species 1
suit_sp1 <- list(suit_sp1)
suit_cube_sp1 <- do.call("c", lapply(suit_sp1, stars::st_as_stars)) %>% setNames(., c("suit"))

# stars cube for Species 2
suit_sp2 <- list(suit_sp2)
suit_cube_sp2 <- do.call("c", lapply(suit_sp2, stars::st_as_stars)) %>% setNames(., c("suit"))

# print(suit_cube_sp1)

## Virtual Data Cube
# bounding box
bb_lux <- st_bbox(tmin)

# from bounding box to polygon
sf_lux <- st_as_sfc(bb_lux) %>% 
  st_sf()

# grid
lux_grid <- st_make_grid(sf_lux, cellsize = .1, n = c(100, 100), what = "polygons", square = FALSE, offset = st_bbox(sf_lux)[c("xmin", "ymin")]) %>% 
  st_as_sf() %>% 
  mutate(id = 1:nrow(.))

# plot raster with grid
# Convert SpatRaster to dataframe
tmin_df <- as.data.frame(tmin$LUX_wc2.1_30s_tmin_1, xy = TRUE, na.rm = TRUE)

ggplot() +
  # Add raster layer
  geom_raster(data = tmin_df, aes(x = x, y = y, fill = LUX_wc2.1_30s_tmin_1)) +
  scale_fill_viridis_c(alpha = 1, begin = 0, end = 1, option = "viridis") +
  # Add grid
  geom_sf(data = lux_grid, color = "black", size = 0.5, fill = NA) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# aggregate by cells, calculating the average of the suitability values within each cell
agg_sp1 <- aggregate(suit_cube_sp1, lux_grid, mean, as_points = TRUE, na.action = na.omit)

agg_sp2 <- aggregate(suit_cube_sp2, lux_grid, mean, as_points = TRUE, na.action = na.omit)

# print(agg_sp1)
cube_sp1sp2 <- c(agg_sp1, agg_sp2) %>%  
  st_redimension() %>% 
  st_set_dimensions(., which = "new_dim", values = c("specie1","specie2"), names = "species") %>%
  st_set_dimensions(., which = "time", values = c("jan","feb"), names = "month") %>% 
  
  setNames(.,"suitability")

print(cube_sp1sp2)

# first, we need to identify the corresponding cell for that location.
which_cell <- st_sf(geometry = st_sfc(st_point(c(5.5, 49.0)), crs = 4326))  %>%  st_join(., lux_grid) 
print(which_cell$id)
# [1] 10

# This object provides, for each species (,,1) and (,,2), the suitability in January and February in cell 10.

# In the first case, it decreased; in the second, it increased.
pull(cube_sp1sp2[,10,,], "suitability")

# We can better visualize this information
values_suit <- pull(cube_sp1sp2[,10,,], "suitability")
# Transform values_suit into a long format without rewriting it
df_long <- melt(values_suit)
colnames(df_long) <- c("cell", "time", "species", "suitability")

# Convert time and species into factors
df_long$time <- factor(df_long$time, labels = c("Jan", "Feb"))
df_long$species <- factor(df_long$species, labels = c("Species 1", "Species 2"))

# plot
ggplot(df_long, aes(x = species, y = suitability, color = species)) +
  geom_point(size = 3) +
  facet_wrap(~time, scales = "free_x", ncol = 2, strip.position = "bottom") + 
  labs(title = "Suitability for species in January and February",
       x = "Month",
       y = "Suitability") +
  scale_color_manual(values = c("Species 1" = "purple2", "Species 2" = "orange")) + 
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),            
    panel.spacing = unit(2, "lines"),          
    axis.line = element_line(color = "black", size = 0.8),
    strip.placement = "outside", 
    legend.position = "right",                 
    panel.border = element_rect(color = "black", fill = NA, size = 0.
    ))
