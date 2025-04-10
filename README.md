# B-Cubed/Suitability Cube

## Description
Scripts to build suitability data cubes

* [Milestone MS20: Code development for predictive habitat suitability modelling](#virtual-suitability-cube-code-development)
* [Milestone MS25: Code testing for predictive habitat suitability](#suitability-cube-code-testing)

### Authors
Rocìo Beatriz Cortès Lobos (Unibo), Matilde Martini (Unibo), Michele Di Musciano (UnivAQ), Duccio Rocchini (Unibo)

### Next steps (Nov. 2025)
* Add temporal dimension ⏲️
* Add an uncertainty measure (```dubicube``` R package)
* Add the Dissimilarity Index to measure the distance between training points and predictors
* Incorporate DeepMaxEnt


# Suitability Cube
Species suitability refers to how favorable an environment is for a species to survive, reproduce, and grow in a specific area and time. It takes into account factors like climate, landscape, and resource availability.


Species Distribution Models (SDMs), also known as Environmental Niche Models (ENMs), are tools that use environmental and species occurrence data to study and predict the distribution of species across time and space. SDMs help identify suitable habitats, forecast the movements of invasive species, and illustrate how species distributions might change due to factors like climate change. They are essential for conservation, allowing us to study how species interact with their environment and make informed decisions to protect biodiversity.

Studying species suitability under different environmental conditions is important for understanding population dynamics, planning conservation actions, and monitoring the effects of climate change and human activities on ecosystems. With this knowledge, we can make better decisions on how to protect habitats and species sustainably.

To facilitate the observation of suitability for multiple species over time and space, we developed a framework that uses **Data Cubes**, multidimensional arrays that organize data in a structured way. In this tutorial, we outline the steps to create a **stars** object, which includes three dimensions: time, space (represented as grid cells), and species, with suitability as the main attribute. Stars objects can be sliced, aggregated along one of the dimensions, and analyzed, making them ideal for studying species suitability.

## Virtual Suitability Cube (Code Development)

For demonstrating the structure of the data cube, we use **virtual species**, which are artificially generated species with known suitability maps based on climate data. The main steps include combining climate data to calculate the suitability of two different species over time and in the same area, then merging these species into a single stars object.

Starting with a time series of climate variables, we combine them to create the suitability for two different virtual species, whose trends over time we want to observe.


<p align="center">
  <img width="610" height="400" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/vsc_page-0001.jpg">
</p>

The two data cubes of species suitability are treated as separate entities, which we can then combine by aggregating them over polygons, creating a vector data cube.

<p align="center">
  <img width="610" height="470" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/vsc_page-0002.jpg">
</p>

This approach makes it easy to visualize and analyze species suitability across time and space for multiple species at once.

### Climatic Data

``` r
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
```


To begin, we will download climate data from [WorldClim](https://www.worldclim.org/).
Climate data are among the datasets from which the suitability of a species can be calculated. 

In this simple example, we will download only 3 variables: minimum temperature, maximum temperature, and precipitation.
The study area is Luxembourg. 

``` r
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
```
<p align="center">
  <img width="500" height="350" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/lux.png">
</p>


Even if there are 12 rasters, each corresponding to a month, the temporal dimension is not expressed in the RasterStack and **needs to be added** by `time` function

``` r
# here, we will simply use numbers from 1 to 12 to represent # the months.
time(tmin) <- time(tmax) <- time(prec) <- 1:12
climate_vars <- c("tmin", "tmax", "prec")

print(tmin)
# class       : SpatRaster 
# dimensions  : 180, 180, 12  (nrow, ncol, nlyr)
# resolution  : 0.008333333, 0.008333333  (x, y)
# extent      : 5.5, 7, 49, 50.5  (xmin, xmax, ymin, ymax)
# coord. ref. : lon/lat WGS 84 (EPSG:4326) 
# source      : LUX_wc2.1_30s_tmin.tif 
# names       : LUX_w~min_1, LUX_w~min_2, LUX_w~min_3, LUX_w~min_4, LUX_w~min_5, LUX_w~min_6, ... 
# min values  :        -3.1,        -3.2,        -0.8,         1.2,         5.6,         8.2, ... 
# max values  :         0.0,        -0.1,         2.6,         4.8,         9.0,        12.1, ... 
# time (raw)  : 1 to 12 
```
<p align="center">
  <img width="600" height="400" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/vsc_pages-to-jpg-0002.jpg">
</p>


### Stars Data Cube
The [stars R package](https://r-spatial.github.io/stars/) provides an infrastructure for managing **data cubes**. Data cubes are multi-dimensional arrays that enable the organization and analysis of large datasets across multiple dimensions, such as time, space, and various environmental variables.

We will create a stars object that contains the climatic variables as attributes and x and y as dimensions. This way, we have everything we need to proceed into a single object.
```r
# first, we create a list that contains all the stacks, to which we apply st_as_stars.
stars_clima <- list(tmin, tmax, prec) %>%
  lapply(st_as_stars) %>%
  do.call("c", .) %>%
  setNames(climate_vars)

print(stars_clima)
# stars object with 3 dimensions and 3 attributes
# attribute(s):
#      Min. 1st Qu. Median      Mean 3rd Qu.  Max.
# tmin  -3.2     0.4    4.4  4.917869     9.4  14.1
# tmax   0.9     5.9   12.5 12.636610    18.9  25.1
# prec  46.0    67.0   75.0 76.737418    84.0 143.0
# dimension(s):
#    from  to offset     delta refsys point x/y
# x       1 180    5.5  0.008333 WGS 84 FALSE [x]
# y       1 180   50.5 -0.008333 WGS 84 FALSE [y]
# time    1  12      1         1     NA    NA  
```

### Virtual Species
The purpose of this brief tutorial is to construct a data structure that allows for the comparison of suitability between two species occupying the same area over time.

As an example, we will consider only two species and two months to observe changes in suitability. The suitability of a real species can be obtained by applying **Species Distribution Models (SDMs)**. But in our case, since we are interested in creating the structure, we will use [virtualspecies](https://borisleroy.com/files/virtualspecies-tutorial.html), which are artificial species randomly generated.

During the creation of these random virtual species, the first step involves generating suitability based on the initial environmental data: in our case, the climatic data.

We will calculate the suitability for each species over the first two months, january and february. 

<p align="center">
  <img width="700" height="450" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/vsc_pages-to-jpg-0003.jpg">
</p>


```r
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
# class      : RasterStack 
# dimensions : 180, 180, 32400, 3  (nrow, ncol, ncell, nlayers)
# resolution : 0.008333333, 0.008333333  (x, y)
# extent     : 5.5, 7, 49, 50.5  (xmin, xmax, ymin, ymax)
# crs        : +proj=longlat +datum=WGS84 +no_defs 
# names      :  tmin,  tmax,  prec 
# min values :  -3.1,   0.9,  53.0 
# max values :   0.0,   4.9, 135.0 

plot(january, col = plasma(500, alpha = 1, begin = 0, end = 1, direction = 1))

```
<p align="center">
  <img width="500" height="350" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/clima_stars.png">
</p>


The following function uses the `generateRandomSp` function from the `virtualspecies` , which creates random species based on the stack of climatic variables and allows us to set some parameters.

Each time the function is executed, a new random species is generated. Therefore, to ensure that we obtain the same species, we need to set a fixed seed using `set.seed()`
```r
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

# same for species 2
suit_sp2_jan <- generate_suitability(january, seed_value = 456)
suit_sp2_feb <- generate_suitability(february, seed_value = 456)
suit_sp2 <- c(suit_sp2_feb, suit_sp2_jan)  

# add time dimension
time(suit_sp2) <- time(suit_sp1) <- 1:2  

# suitability for Species 1
plot(suit_sp1, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
title("Suitability of Species 1 in Jan-Feb", outer=TRUE, line=-0.9)
```

<p align="center">
  <img width="500" height="350" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/suit_sp1.png">
</p>

```r
# suitability for Species 2
plot(suit_sp2, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
title("Suitability of Species 2 in Jan-Feb", outer=TRUE, line=-0.9)
```

<p align="center">
  <img width="500" height="350" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/suit_sp2.png">
</p>

For each species, we will create a data cube: the dimensions remain x, y, and time, but this time the attribute is the suitability only, with values ranging from 0 to 1, which characterizes each species.
```r
# stars cube for Species 1
suit_sp1 <- list(suit_sp1)
suit_cube_sp1 <- do.call("c", lapply(suit_sp1, stars::st_as_stars)) %>% setNames(., c("suit"))

# stars cube for Species 2
suit_sp2 <- list(suit_sp2)
suit_cube_sp2 <- do.call("c", lapply(suit_sp2, stars::st_as_stars)) %>% setNames(., c("suit"))
```
```r
print(suit_cube_sp1)
# stars object with 3 dimensions and 1 attribute
# attribute(s):
#      Min.   1st Qu.    Median      Mean   3rd Qu. Max.
# suit     0 0.1493433 0.5079458 0.4706163 0.7610167    1
# dimension(s):
#      from  to offset     delta refsys x/y
# x       1 180    5.5  0.008333 WGS 84 [x]
# y       1 180   50.5 -0.008333 WGS 84 [y]
# time    1   2      1         1     NA 

print(suit_cube_sp2)
# stars object with 3 dimensions and 1 attribute
# attribute(s):
#       Min.     1st Qu.     Median      Mean   3rd Qu. Max.
# suit     0 0.009498931 0.06620133 0.1808423 0.3206807    1
# dimension(s):
#      from  to offset     delta refsys x/y
# x       1 180    5.5  0.008333 WGS 84 [x]
# y       1 180   50.5 -0.008333 WGS 84 [y]
# time    1   2      1         1     NA   
```

### Virtual Data Cube

We want to incorporate both species into the same object.

To combine the two species into a single object, we need to collapse the x and y dimensions into one, creating a grid that represents individual cells within the study area. This grid allows us to reduce the dimensions of the stars object and facilitates the transition to a **vector data cube**.

In this process, pixels in the raster are grouped based on their spatial intersection with a set of vector geometries. Each group is then reduced to a single value using an aggregation function, such as the mean or maximum. The result is a one-dimensional sequence of feature geometries defined in space.

By creating the grid, it is possible to reduce two dimensions (x and y) into a single one (cell): this way, we can observe the suitability over time for both species within the defined polygons. As a result, "species" becomes an additional dimension in our data structure.
``` r
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
```

<p align="center">
  <img width="350" height="400" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/grid_lux_coordinates.png">
</p>

```r
# aggregate by cells, calculating the average of the suitability values within each cell
agg_sp1 <- aggregate(suit_cube_sp1, lux_grid, mean, as_points = TRUE, na.action = na.omit)

agg_sp2 <- aggregate(suit_cube_sp2, lux_grid, mean, as_points = TRUE, na.action = na.omit)

print(agg_sp1)
# stars object with 2 dimensions and 1 attribute
# attribute(s):
#               Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
# suit  0.0009594404 0.1817593 0.5260862 0.4765746 0.7482857 0.9380438   36
# dimension(s):
#      from  to offset delta refsys point                                   values
# x       1 297     NA    NA WGS 84 FALSE            POLYGON ((5.45 49.02887, ...
# time    1   2      1     1     NA   NA                                   NULL

print(agg_sp2)
# stars object with 2 dimensions and 1 attribute
# attribute(s):
#               Min.    1st Qu.     Median      Mean   3rd Qu.      Max. NA's
# suit  2.707365e-05 0.01614575 0.07359824 0.1784444 0.3089694 0.9104424   36
# dimension(s):
#     from  to offset delta refsys point                                 values
# x       1 297     NA    NA WGS 84 FALSE POLYGON ((5.45 49.02887, ...
# time    1   2      1     1     NA    NA                                  NULL
```
Finally, we merge the two cubes. Now the dimensions are 3 again: cell, time and species. The attribute represents the species' suitability over time and space.

``` r
cube_sp1sp2 <- c(agg_sp1, agg_sp2) %>%  
  st_redimension() %>% 
  st_set_dimensions(., which = "new_dim", values = c("specie1","specie2"), names = "species") %>%
  st_set_dimensions(., which = "time", values = c("jan","feb"), names = "month") %>% 
  
  setNames(.,"suitability")

print(cube_sp1sp2)
# stars object with 3 dimensions and 1 attribute
# attribute(s):
#                      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
# suitability  2.707365e-05 0.0440898 0.2516058 0.3275095 0.5842543 0.9380438   72
# dimension(s):
#        from  to refsys point                                                        values
# x          1 297 WGS 84 FALSE POLYGON ((5.45 49.02887, ...,...,POLYGON ((7.05 50.41451, ...
# month      1   2     NA    NA                                                      jan, feb
# species    1   2     NA    NA                                              specie1, specie2

```

Let's suppose we have a specific point in space: we want to investigate the suitability of the two species at that location
With `pull` function you can extract the attributes in a specific point of space and time.
```r
# first, we need to identify the corresponding cell for that location.
which_cell <- st_sf(geometry = st_sfc(st_point(c(5.5, 49.0)), crs = 4326))  %>%  st_join(., lux_grid) 
print(which_cell$id)
# [1] 10

# this object provides, for each species (,,1) and (,,2), the suitability in January and February in cell 10
# in the first case, it decreased; in the second, it increased.
pull(cube_sp1sp2[,10,,], "suitability")
# , , 1
#          [,1]      [,2]
# [1,] 0.9036776 0.8826196
#
# , , 2
#            [,1]        [,2]
# [1,] 0.001298286 0.007460382

## we can better visualize this information this way:
values_suit <- pull(cube_sp1sp2[,10,,], "suitability")
# transform values_suit into a long format without rewriting it
df_long <- melt(values_suit)
colnames(df_long) <- c("cell", "time", "species", "suitability")

# convert time and species into factors
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
    panel.border = element_rect(color = "black", fill = NA, size = 0.))
```
<p align="center">
  <img width="550" height="330" src="https://github.com/rociobeatrizc/virtual_suitability_cube/blob/main/images/boxplot.png">
</p>


The following is an example of suitability over the course of a whole year (12 months) for 5 different species in two distinct locations in the Luxembourg area.

<p align="center">
  <img width="950" height="550" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/2cells.jpg">
</p>

We can imagine the structure of the cube in this way, where the temporal layers are 12

<p align="center">
  <img width="450" height="450" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/5sp_cube.jpg">
</p>

## Suitability Cube (Code testing)

Code testing with real data from Dutch Vegetation Database. 
Here, we want to test the stars data cube structure for multiple real species. The goal is to predict the suitability of those species in a different area from the training one. 

Steps are: 
* Collecting climatic data of The Netherlands (one month) and creating a stars data cube with them
* Collecting only-presence data for 5 species in Netherlands from 2000 to 2017 and splitting them by species name
* Using the environmental information and the occurrences, training a MaxEnt model for each species
* Collecting climatic data of Belgium (one month) and creating a stars data cube with them
* Predicting the suitability of the 5 species in the new area using models trained before
* Aggregating the suitability maps over polygons to synthetize all the informations into a single data cube

```r
# imports
library(ggplot2)
library(terra)
library(sf)
library(raster)
library(viridis)
library(tidyverse)
library(stars)
library(geodata)
library(purrr)
library(reshape2)
library(dplyr)
library(enmSdmX)
```
### Climatic data for the training area (The Netherlands)

Data collection and aggregation in a data cube that contains the climatic variables as attributes and x, y and time. 

``` r
## Download climatic data for the Netherlands
# worldclim_country function
# The res paramater is the resolution, with valid values as 10,5,2.5 and 0.5 
# The var parameter is the variable name, with valid values as tmin, tmax, tavg, prec, wind, vapr and bio.
# The worldclim_country function will return 12 layers of raster (in a SpatRaster), for each variable. 
tmin <- worldclim_country("Nld", "tmin", path=tempdir(), res = 0.5, version = "2.1")
tmax <- worldclim_country("Nld", "tmax", path=tempdir(), res = 0.5, version = "2.1")
prec <- worldclim_country("Nld", "prec", path=tempdir(), res = 0.5, version = "2.1")
tavg <- worldclim_country("Nld", "tavg", path=tempdir(), res = 0.5, version = "2.1")
wind <- worldclim_country("Nld", "wind", path=tempdir(), res = 0.5, version = "2.1")

# add time
time(tmin) <- time(tmax) <- time(prec) <- time(tavg) <- time(wind) <- 1:12

climate_vars <- c("tmin", "tmax", "prec", "tavg", "wind")

# first, we create a list that contains all the stacks, to which we apply st_as_stars.
stars_clima <- list(tmin, tmax, prec, tavg, wind) %>%
  lapply(st_as_stars) %>%
  do.call("c", .) %>%
  setNames(climate_vars)

print(stars_clima)

# stars object with 3 dimensions and 5 attributes
# attribute(s), summary of first 1e+05 cells:
#  Min. 1st Qu. Median         Mean 3rd Qu. Max.  NA's
# tmin  -0.8    -0.4   -0.1  0.007868939     0.2  2.6 63985
# tmax   3.7     4.3    4.6  4.718736640     5.0  6.4 63985
# prec  58.0    66.0   68.0 67.980036096    70.0 80.0 63985
# tavg   1.6     1.9    2.2  2.363004291     2.6  4.5 63985
# wind   4.0     4.6    5.2  5.363057057     6.0  7.8 63985
#
# dimension(s):
#     from  to offset     delta refsys point x/y
# x       1 540      3  0.008333 WGS 84 FALSE [x]
# y       1 420     54 -0.008333 WGS 84 FALSE [y]
# time    1  12      1         1     NA    NA    

# training month
clima_may <- stars_clima %>% slice("time", 5) 

# predictors
clima_train <- rast(clima_may) %>% 
  setNames(climate_vars)

plot(clima_train)
```

<p align="center">
  <img width="550" height="450" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/netherlands_clima.png">
</p>


### Occurrences

From Dutch Vegetation Database: [GBIF Occurrence Download](https://doi.org/10.15468/dl.hynkda)
Occurrences of 5 species, collected from 2000 to 2017
* Galium verum L.
* Ophrys apifera Huds.
* Paris quadrifolia L.
* Chrysosplenium alternifolium L.
* Anemone nemorosa L.

They will be use for training MaxEnt model
``` r
## occurrences

# upload the dataset
# read txt
occ <- read.delim("occurrence.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  select(scientificName, decimalLatitude, decimalLongitude, year)

## plot
# Filter years from 2000 to 2017
occ_counts <- occ %>%
  filter(year >= 2000 & year <= 2017) %>%
  group_by(scientificName, year) %>%
  summarise(count = n(), .groups = "drop")

# Heatmap to see the occurrences of each species per year
ggplot(occ_counts, aes(x = year, y = scientificName, fill = count)) +
  geom_tile() +
  scale_fill_gradientn(colors = mako(10, direction = -1), 
                       breaks = seq(0, max(occ_counts$count, na.rm = TRUE), by = 150), 
                       labels = seq(0, max(occ_counts$count, na.rm = TRUE), by = 150)) +
  labs(title = "occurrences per species per year",
       x = "year",
       y = "species",
       fill = "occurrences") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed(ratio = 1.5)
```

<p align="center">
  <img width="800" height="550" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/heatmap_page-0001.jpg">
</p>

``` r
# filter year
occ <- occ %>% 
  filter(year >= 2000 & year <= 2017) %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

head(occ)
# scientificName decimalLatitude decimalLongitude year
# 1     Galium verum L.        51.82729          3.95035 2017
# 2     Galium verum L.        52.87226          4.71202 2017
# 3     Galium verum L.        50.84668          5.77765 2015
# 4     Galium verum L.        52.53752          6.47736 2016
# 5 Anemone nemorosa L.        52.01382          6.07912 2013
# 6 Anemone nemorosa L.        50.88266          5.82037 2014

# split dataset by species with the function split_species_data
species <- split_species_data(occ)

typeof(species)
# [1] "list"

names(species)
# [1] "Anemone nemorosa L."             "Chrysosplenium alternifolium L." "Galium verum L."                
# [4] "Ophrys apifera Huds."            "Paris quadrifolia L."  
```
### MaxEnt model for each species
The function takes predictors and occurrences to build a MaxEnt model (enmSdmX R package) for each species.
Models are then saved and can be used for predicting in a new area. The output contains also prediction maps. 

``` r
# creating SDMs with MaxEnt for many species in the same area, with the same predictors
sdms <- create_sdm_for_species_list(species, clima_train, background_points = 10000, predictors = names(clima_train))

# plot suitability map for Anemone nemorosa L. in The Netherlands
plot(sdms$predictions$`Anemone nemorosa L.`)
```

<p align="center">
  <img width="500" height="450" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/anemone.png">
</p>


### New area: Belgium
Based on the models previously trained, let's check the suitability of our species in a different area, i.e. Belgium.

```r
# predictions for the same species but in another area (Belgium)
# same climatic variables
tmin_b <- worldclim_country("Bel", "tmin", path=tempdir(), res = 0.5, version = "2.1")
tmax_b <- worldclim_country("Bel", "tmax", path=tempdir(), res = 0.5, version = "2.1")
prec_b <- worldclim_country("Bel", "prec", path=tempdir(), res = 0.5, version = "2.1")
tavg_b <- worldclim_country("Bel", "tavg", path=tempdir(), res = 0.5, version = "2.1")
wind_b <- worldclim_country("Bel", "wind", path=tempdir(), res = 0.5, version = "2.1")

# time
time(tmin_b) <- time(tmax_b) <- time(prec_b) <- time(tavg_b) <- time(wind_b) <- 1:12

# stars object of Belgium
stars_clima_bel <- list(tmin_b, tmax_b, prec_b, tavg_b, wind_b) %>%
  lapply(st_as_stars) %>%
  do.call("c", .) %>%
  setNames(climate_vars)

# may for Belgium (we want to see suitability predictions for the same month)
clima_may_bel <- stars_clima_bel %>% slice("time", 5) 

# predictors in may
clima_train_bel_may <- rast(clima_may_bel) %>% 
  setNames(climate_vars)

# predicting suitability in a new area for the same species using the models trained previously 
new_predictions_may <- predict_sdm_for_new_area(sdms$models, clima_train_bel_may)

# the output is a data cube with suitability as attribute
print(new_predictions_may)

# stars object with 3 dimensions and 1 attribute
# attribute(s):
#              Min.    1st Qu.    Median      Mean   3rd Qu.      Max.  NA's
# suit  5.867529e-13 0.03843691 0.1433858 0.2354462 0.3725311 0.9999935 67380
# dimension(s):
#        from  to offset     delta refsys                                       values x/y
# x          1 480    2.5  0.008333 WGS 84                                         NULL [x]
# y          1 360     52 -0.008333 WGS 84                                         NULL [y]
# species    1   5     NA        NA     NA Anemone nemorosa L.,...,Paris quadrifolia L. 

```
### Aggregation
The last step is the polygon-based aggregation. 
In line with the structure proposed by the B-Cubed framework, we aim to summarize species suitability information within a grid covering the study area.

```r
## aggregation steps
# bounding box
bbox <- st_bbox(tmin_b)
sf_bel <- st_as_sfc(bbox) %>% 
  st_sf()

# grid
bel_grid <- st_make_grid(sf_bel, cellsize = .1, n = c(50, 50), what = "polygons", square = FALSE, offset = st_bbox(sf_bel)[c("xmin", "ymin")]) %>% 
  st_as_sf() %>% 
  mutate(id = 1:nrow(.))

# plot raster with grid
# convert SpatRaster to dataframe
tmin_df <- as.data.frame(tmin_b$BEL_wc2.1_30s_tmin_1, xy = TRUE, na.rm = TRUE)
tmin_df
p <- ggplot() +
  geom_raster(data = tmin_df, aes(x = x, y = y, fill = BEL_wc2.1_30s_tmin_1)) +
  scale_fill_viridis_c(alpha = 1, begin = 0, end = 1, option = "viridis") +
  geom_sf(data = bel_grid, color = "black", size = 0.5, fill = NA, alpha = 0.2) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )
plot(p)
```

<p align="center">
  <img width="430" height="350" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/grid_bel.png">
</p>


```r
# aggregating suitability values of the species over polygons for a given area
stars_predictions_aggregated <- aggregate_suitability(new_predictions_may, bel_grid)

# output: data cube
print(stars_predictions_aggregated)

# stars object with 2 dimensions and 1 attribute
# attribute(s):
#                     Min.    1st Qu.    Median     Mean  3rd Qu.      Max. NA's
# suitability  1.540887e-07 0.05111451 0.1581302 0.234773 0.370952 0.9726111  930
# dimension(s):
#        from   to refsys point                                                        values
# x          1 1494 WGS 84 FALSE POLYGON ((2.45 49.02887, ...,...,POLYGON ((6.55 51.97335, ...
# species    1    5     NA    NA                  Anemone nemorosa L.,...,Paris quadrifolia L


# first, we need to identify the corresponding cell for that location.
which_cell <- st_sf(geometry = st_sfc(st_point(c(4.3517, 50.8503)), crs = 4326))  %>%  st_join(., bel_grid) 
print(which_cell$id)
# [1] 695

# extract suitability values
pull(stars_predictions_aggregated[,695,], "suitability")
#           [,1]       [,2]       [,3]      [,4]       [,5]
# [1,] 0.00129283 0.03917683 0.02307574 0.3141483 0.07493015
```


## References
* virtualspecies, an R package to generate virtual species distributions (2018), Leroy B. et al., [DOI](https://doi.org/10.1111/ecog.01388)
* Meyer H, Milà C, Ludwig M, Linnenbrink J, Schumacher F (2024). CAST: 'caret' Applications for Spatial-Temporal Models. R package version 1.0.2, [https://hannameyer.github.io/CAST/](https://github.com/HannaMeyer/CAST)
* Gilardi A, Lovelace R (2024). osmextract: Download and Import Open Street Map Data Extracts. R package version 0.5.1.900, [https://github.com/ropensci/osmextract](https://docs.ropensci.org/osmextract/).
* A practical overview of transferability in species distribution modeling, Werkowska W., Marquez A., Real R., Acevedo P., 2017,  [DOI](https://cdnsciencepub.com/doi/abs/10.1139/er-2016-0045)
* Effect of roadside bias on the accuracy of predictive maps produced by bioclimatic models, Kadmon R., Farber O., Danin A., 2004, [DOI](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/02-5364)
* The virtual ecologist approach: simulating data and observers, Zurell et al., 2010, [DOI](https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0706.2009.18284.x)
* The n-dimensional hypervolume, Blonder et al., 2014, [DOI](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12146)
* The cumulative niche approach: a framework to assess the performance of ecological niche model projections, Arlè et al., 2024, [DOI](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.11060)
*  Smith, A.B., Murphy, S.J., Henderson, D., and Erickson, K.D. 2023. Including imprecisely georeferenced specimens improves accuracy of species distribution models and estimates of niche breadth. Global Ecology and Biogeography In press
* Pebesma E, Bivand R (2023). Spatial Data Science: With applications in R. Chapman and Hall/CRC, London. doi:10.1201/9780429459016, https://r-spatial.org/book/.
*	Meyer, Hanna, and Edzer Pebesma. "Predicting into unknown space? Estimating the area of applicability of spatial prediction models." Methods in Ecology and Evolution 12.9 (2021): 1620-1633.
*	Steven J. Phillips, Robert P. Anderson, Robert E. Schapire, Maximum entropy modeling of species geographic distributions, Ecological Modelling, Volume 190, Issues 3–4, 2006, Pages 231-259, ISSN 0304-3800, https://doi.org/10.1016/j.ecolmodel.2005.03.026.
*	Guisan, Antoine, and Wilfried Thuiller. "Predicting species distribution: offering more than simple habitat models." Ecology letters 8.9 (2005): 993-1009.
