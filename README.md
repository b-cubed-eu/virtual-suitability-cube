# B-Cubed/Suitability Cube

Milestone 4



# Virtual Suitability Cube
Species suitability refers to how favorable an environment is for a species to survive, reproduce, and grow in a specific area and time. It takes into account factors like climate, landscape, and resource availability.


Species Distribution Models (SDMs), also known as Environmental Niche Models (ENMs), are tools that use environmental and species occurrence data to study and predict the distribution of species across time and space. SDMs help identify suitable habitats, forecast the movements of invasive species, and illustrate how species distributions might change due to factors like climate change. They are essential for conservation, allowing us to study how species interact with their environment and make informed decisions to protect biodiversity.

Studying species suitability under different environmental conditions is important for understanding population dynamics, planning conservation actions, and monitoring the effects of climate change and human activities on ecosystems. With this knowledge, we can make better decisions on how to protect habitats and species sustainably.

To facilitate the observation of suitability for multiple species over time and space, we developed a framework that uses **Data Cubes**, multidimensional arrays that organize data in a structured way. In this tutorial, we outline the steps to create a **stars** object, which includes three dimensions: time, space (represented as grid cells), and species, with suitability as the main attribute. Stars objects can be sliced, aggregated along one of the dimensions, and analyzed, making them ideal for studying species suitability.

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

## Climatic Data

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
plot(tmin)
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

## Stars Data Cube
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

## Virtual Species
The purpose of this brief tutorial is to construct a data structure that allows for the comparison of suitability between two species occupying the same area over time.

As an example, we will consider only two species and two months to observe changes in suitability. The suitability of a real species can be obtained by applying **Species Distribution Models (SDMs)**. But in our case, since we are interested in creating the structure, we will use [virtualspecies](https://borisleroy.com/files/virtualspecies-tutorial.html), which are artificial species randomly generated.

During the creation of these random virtual species, the first step involves generating suitability based on the initial environmental data: in our case, the climatic data.

We will calculate the suitability for each species over the first two months, january and february. 


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
  set.seed(seed_value)  # Assicura che la specie virtuale sia riproducibile
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
## Virtual Data Cube

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

```r
# first, we need to identify the corresponding cell for that location.
which_cell <- st_sf(geometry = st_sfc(st_point(c(5.5, 49.0)), crs = 4326))  %>%  st_join(., lux_grid) 
print(which_cell$id)
# [1] 10

# This object provides, for each species (,,1) and (,,2), the suitability in January and February in cell 10.

# In the first case, it decreased; in the second, it increased.
pull(cube_sp1sp2[,10,,], "suitability")
# , , 1
#          [,1]      [,2]
# [1,] 0.9036776 0.8826196
#
# , , 2
#            [,1]        [,2]
# [1,] 0.001298286 0.007460382

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


## References
* virtualspecies, an R package to generate virtual species distributions (2018), Leroy B. et al., [DOI](https://doi.org/10.1111/ecog.01388)
* Pebesma E, Bivand R (2023). Spatial Data Science: With applications in R. Chapman and Hall/CRC, London. DOI:10.1201/9780429459016, [https://r-spatial.org/book/](https://r-spatial.org/book/)
* Meyer H, Milà C, Ludwig M, Linnenbrink J, Schumacher F (2024). CAST: 'caret' Applications for Spatial-Temporal Models. R package version 1.0.2, [https://hannameyer.github.io/CAST/](https://github.com/HannaMeyer/CAST)
