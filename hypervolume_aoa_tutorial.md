``` r
# load packages
library(CAST)
library(caret)
library(sf)
library(ClimDatDownloadR)
library(devtools)
library(raster)
library(viridis)
library(virtualspecies)
library(ggplot2)
library(tidyverse)
library(terra)
library(ade4)
library(hypervolume)
library(gdata)
library(ggmap)
library(osmdata)
library(osmextract)

# set working directory
setwd("my/path")
```
## Input data
Carichiamo i file vettoriali: l'area di studio e il reticolo stradale da Open Street Map ([osmextract R package](https://docs.ropensci.org/osmextract/))

``` r
# upload shapefile
aoi <- st_read("aoi.shp")

# bounding box 
aoi_bb <- st_bbox(aoi)

# from Open Street Map select type of roads: primary, secondary, tertiary (paths)
ht_secondary <- "secondary"

# download roads from OSM: our region is Abruzzo, Italy
osm_aoi <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_aoi_roads <- osm_aoi[osm_aoi$highway %in% ht_secondary, ]
plot(osm_aoi_roads$geometry)

```
The input data for virtual species are environmental spatial data (raster data).

**CHELSA** (Climatologies at high resolution for the earth’s land surface areas) is a very high resolution (30 arc sec, ~1km) global downscaled climate data set: it is built to provide free access to high resolution climate data for research and application, and is constantly updated and refined.

[Bioclimatic variables](https://chelsa-climate.org/bioclim) are derived variables developed for species distribution modeling and related ecological applications.

You can directly download CHELSA data into R.

``` r

# download bioclimatic variables from CHELSA
Chelsa.Clim.download(
  # Starting from the workind directory, specify the path
  #  save.location = "my/path",
  
  # 'bio' contains all the bioclimatic variables
  parameter = "bio",
  
  # some variables are chosen from the 19 available
  bio.var = c(1, 7, 13, 14),
  
  # version
  version.var = "2.1",
  
  # cropping along the area of interest
  clipping = TRUE,
  clip.shapefile = aoi,
  
  # insert the coordinates of the area of interest (bounding box)
  clip.extent = aoi_bb,
  
  # buffer, if needed
  # buffer = 3,
  
  # other commands
  convert.files.to.asc = FALSE,
  stacking.data = TRUE,
  combine.raw.zip = FALSE,
  delete.raw.data = FALSE,
  save.bib.file = TRUE
)


##  upload bioclimatic variables
# string containing the names of raster files
# NB: The following pattern means that we want to read *.tif files which are
# included in a directory that contains the word 'clipped'
# (that's the style adopted by ClimDatDownloadR package)
# Feel free to modify it if you have different requirements
rastlist <- list.files(path ="my/path/bio/ChelsaV2.1Climatologies/clipped_2024-09-18_10-46-34", pattern = "CHELSA", full.names = TRUE)

# using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# change data names
names(mydata) <- c("bio01", "bio07", "bio13", "bio14")

# crop and mask by region borders
aoi_sp <- sf::as_Spatial(aoi)
mydata <- mydata %>% crop(., aoi_sp) %>% mask(., aoi_sp)

# original data: will be useful later
mydata_backup <- mydata
```
## Virtual Species
Generating random species from known environmental data allows controlling the factors that can influence the distribution of real data. Questo può essere fatto con il pacchetto virtualspecies.  
To create a series of occurrence points for a species, it is necessary to go through 3 steps. All'interno di ciascuno step si possono regolare dei parametri il cui significato è spiegato nel tutorial virtualspecies:

1. By intersecting bioclimatic data, the first output obtained is a **suitability map**
2. The suitability map is converted into a binary **presence/absence map** through a probability function (logistic curve) that associates the suitability value with the probability of finding the virtual species for each pixel. This subset of the environmental niche that is actually occupied by the species corresponds to the realized niche (Hutchinson, 1957).
3. The third step consists of generating, within the presence/absence raster, a series of occurrence points. 


``` r
## Random Virtual Species: run every time you want to create a virtual species.

## step 1: suitability map generation
random.sp <- generateRandomSp(raster.stack = mydata,
                              convert.to.PA = FALSE,
                              # how to combine response functions
                              species.type = "multiplicative",
                              # random approach between PCA and response function
                              approach = "response",
                              # response function
                              relations = "gaussian",
                              # realistic species
                              realistic.sp = TRUE,
                              plot = FALSE)

## step 2: Presence/Absence: requires defining the parameters alpha, beta, and species prevalence
new.pres <-convertToPA(random.sp,
                       beta = "random",
                       alpha = -0.05, plot = FALSE,
                       species.prevalence = 0.01)

## step 3: occurences
presence.points <- sampleOccurrences(new.pres,
                                     n = 200,
                                     type = "presence only",
                                     sample.prevalence = 0.9,
                                     error.probability = 0,
                                     detection.probability = 1,
                                     correct.by.suitability = TRUE,
                                     plot = FALSE)
```

<p align="center">
  <img width="650" height="250" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/vsc_page-0003.jpg">
</p>


## Hypervolume
Hutchinson defined an ecological niche as an n-dimensional volume in
the environmental space where a species can maintain a viable population and persist along
time (Sillero et al., 2021). The Hutchinsonian niche, or n-dimensional environmental space, is
defined as the hypervolume (Blonder et al., 2014). 

The calculation of the hypervolume is performed using the 'hypervolume_gaussian' function from the R package hypervolume, developed by Blonder et al. `hypervolume_gaussian` constructs a hypervolume by building a gaussian kernel density estimate on an adaptive grid of random points wrapped around the original data points. 

However, since we aim to calculate not just one but multiple hypervolumes as we increase and accumulate occurrences, `hypervolume_gaussian` will be at the core of a function that, starting from a single random occurrence, progressively increases the number of occurrences (and their corresponding bioclimatic variables) and calculates the hypervolume. This increment, a subsample of the original sample, is random and without replacement.

```r
## preliminary Steps for Niche Analysis

# Z transform for hypervolume building
for (i in 1:nlayers(mydata)){
  mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd') 
}

# the raster of occurrences is transformed into a dataset, from which the rows satisfying both conditions Real = 1 and Observed = 1 are preserved
raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]

# the environmental variables are associated with the occurrences using their coordinates
values_occ <- mydata %>% rasterToPoints() %>% as.data.frame()
filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))
occurrences_values <- filtered_occ[,-c(1:2, 7:8)]

## Functions for Hypervolume

# Hypervolume: just the hypervolume value from hypervolume_gaussian function
hyp_calc <- function(data) {
  hv_occ <- hypervolume_gaussian(data)
  return(hv_occ@Volume)
}

# Function to build the accumulation curve with random increment in occurrences
acc_curve <- function(x, no) {
  # Starts with a random row
  fx <- x %>% 
    sample_n(size = 1) 
  
  ipervolumi <- 0
  num_occurrences <- 0
  
  for (i in 1:1000) {
    
    # To the initial value (a row)
    # Random values are selected
    # They are bound to fx
    # Unique values are kept
    fx <- x %>% 
      sample_n(size = no) %>% 
      bind_rows(fx) %>% 
      distinct()
    
    # Hypervolume per subset
    hv <- hyp_calc(fx)
    
    # Save hypervolume & number of occurrences
    ipervolumi <- c(ipervolumi, hv)
    num_occurrences <- c(num_occurrences, nrow(fx))
    
    # Condition
    # Stop when the subset has the same number of occurrences as the original set
    if(nrow(fx) == nrow(x)) {
      break
    }
  }
  
  result <- bind_cols(iperv = ipervolumi, n_occ = num_occurrences)
  return(list(result))
}
```
## Roadside Bias 
To obtain a subsample that shows sampling bias, we first need to associate each point in the original dataset with the distance to the nearest road. Instead of calculating this distance for each point, we chose to rasterize the distances for simplicity: this way, each pixel will have a value representing the distance to the nearest road. Points are then projected onto this raster, and their respective distance values are extracted 

``` r

## Roadside bias
# create raster with distances from roads
roads_vect <- terra::vect(osm_aoi_roads$geometry)

# turn into SpatRaster object
raster_roads <- as(mydata_backup[[1]], "SpatRaster")

# rasterize distances
r <- terra::rasterize(roads_vect, raster_roads)
d <- distance(r, unit = "km") 

## plot: distance from roads
d_rast <- d %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
raster_df_dist <- as.data.frame(rasterToPoints(d_rast), xy = TRUE)
value_column <- names(raster_df_dist)[3]

ggplot() +
  # Add raster
  geom_raster(data = raster_df_dist, aes_string(x = "x", y = "y", fill = value_column)) +
  scale_fill_viridis_c(alpha = 1, begin = 0, end = 1) +  
  # Add roads
  geom_sf(data = osm_aoi_roads$geometry, color = "black", size = 0.5) +
  theme_bw() +
  theme_minimal() +
  labs(title = "Distance from Roads",
       fill = "Distance (km)") +
  coord_sf() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
```
<p align="center">
  <img width="420" height="450" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/distance_poster.png">
</p>


As the distance from the road network increases, the probability of sampling a species will decrease. We can simulate this situation with the following probability function

$$P(sampling) = 1 - \dfrac{log(c\cdot min.dist)}{log(c \cdot max(min.dist))}$$

We used this function to transform distances into probability values: for each distance value, there is a probability for a point to be sampled. Once again, we obtain a raster, making it easy to extract the probability value.

After setting a probability threshold, we selected the points that fall on pixels containing values equal to or above the threshold. We decided to select the points that have a 100% probability of being sampled, i.e., those located on roads or in their immediate surroundings

``` r

## extract distances
d_raster <- d %>% raster()
distances <- d_raster %>%  as.data.frame()

## sampling probability function: simulation of the lazy sampler
c <- 1
sampling_prob <- 1-(((log(c*distances))/(log(max(c*distances)))))
sampling_prob <- as.data.frame(sampling_prob)

# some values are: Inf. Replace those values with 1
sampling_prob[sampling_prob == Inf] <- 1
sampling_prob[sampling_prob > 1] <- 1

# new raster with probability to be sampled instead of distances
prob_raster <- classify(d, cbind(values(d), sampling_prob))

## plot: sampling probability
prob_r <- prob_raster %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
raster_df_prob <- as.data.frame(rasterToPoints(prob_r), xy = TRUE)
value_column <- names(raster_df_prob)[3]

ggplot() +
  # Add raster
  geom_raster(data = raster_df_prob, aes_string(x = "x", y = "y", fill = value_column)) +
  scale_fill_viridis_c(alpha = 1, begin = 0, end = 1) +  
  # Add roads
  geom_sf(data = osm_aoi_roads$geometry, color = "black", size = 0.5) +
  theme_bw() +
  theme_minimal() +
  labs(title = "Sampling Probability",
       fill = "Probability") +
  coord_sf() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

```
<p align="center">
  <img width="420" height="450" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/sampling_prob_poster.png">
</p>


``` r

## Occurrences as points
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")


# Add probability value
points_biased <- coord_occ %>%
  cbind(terra::extract(prob_raster, ., ID = FALSE)) %>%
  subset(.$layer == 1)


## Hypervolume of occurrences (random sampled: null model) 
# Num. simulations each species
num_sim <- 3

# Set stop point according to the number of biased occurrences: same sampling effort
nrow(points_biased)
nrow(occurrences_values)
stop <-  ceiling(nrow(points_biased) + 0.2 * (nrow(points_biased)))

# Random subsample of occurrences from null model: 20%
occurrences_values <- occurrences_values[sample(nrow(occurrences_values), stop), ]

## Plot: map with unbiased-biased points
# Index
indices <- rownames(occurrences_values)
indices <- as.numeric(indices)
filtered_coord_occ <- coord_occ[indices, ]

par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 

plot(prob_r, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
points(filtered_coord_occ, cex = 0.6)
points(points_biased, col = "red", cex = 0.6)
legend("topright", legend = c("Unbiased", "Biased"), col = c("black", "red"), pch = 19, cex = 0.8,
       xpd = TRUE, y.intersp = 0.8)

dev.off()

# List with the occurrences we want to test
hyp_steps <- c(seq(from = 40, to = stop, by = 30), stop)

# Empty list 
all_sim <- list()

# For cycle for simulations
for (sim in 1:num_sim) {
  
  list_output <- list()
  
  for (i in seq_along(hyp_steps)) {
    d_hyp <- acc_curve(occurrences_values, hyp_steps[i])
    list_output[[i]] <- d_hyp[[1]]
    }
  
  all_sim[[sim]] <- list_output
  
}

# All simulations in one df
combined_df <- do.call(rbind, lapply(seq_along(all_sim), function(sim) {
  do.call(rbind, lapply(all_sim[[sim]], function(df) {
    df$sim <- sim
    df
  }))
}))


# Mean predictions (LOESS): x sequence 
x_seq <- seq(min(combined_df$n_occ), max(combined_df$n_occ), length.out = 100)

# Mean predictions: LOESS method
loess_predictions <- lapply(unique(combined_df$n_occ), function(n) {
  preds <- sapply(all_sim, function(lista) {
    loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
    predict(loess_fit, newdata = data.frame(n_occ = n))
  })
  
  data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
  
})


# Mean df
pred_mean <- do.call(rbind, loess_predictions)

# Plot: unbiased hypervolume
ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = pred_mean, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  labs(title = "Mean (unbiased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()

## Hypervolume of biased occurrences (road driven: biased sampling)
biased_df <- points_biased %>%
  as.data.frame() %>%
  .[,-c(5:8)]


# Stop
stop_biased <- nrow(biased_df)
hyp_steps_b <- c(seq(from = 20, to = stop_biased, by = 20), stop_biased)

# Empty list
all_sim_b <- list()

# For cycle for simulations
for (sim in 1:num_sim) {
  
  list_output_b <- list()
  
  for (i in seq_along(hyp_steps_b)) {
    d_hyp <- acc_curve(biased_df, hyp_steps_b[i])
    list_output_b[[i]] <- d_hyp[[1]]
    }
  
  all_sim_b[[sim]] <- list_output_b
  
}


# Combined df
combined_df_biased <- do.call(rbind, lapply(seq_along(all_sim_b), function(sim) {
  
  do.call(rbind, lapply(all_sim_b[[sim]], function(df) {
    df$sim <- sim
    df
    
  }))
  
}))


# Mean LOESS
loess_predictions_biased <- lapply(unique(combined_df_biased$n_occ), function(n) {
  
  preds <- sapply(all_sim_b, function(lista) {
    loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
    predict(loess_fit, newdata = data.frame(n_occ = n))
  })
  
  data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
  
})


# Mean in one df
pred_mean_b <- do.call(rbind, loess_predictions_biased)

## Plot: biased hypervolume
ggplot() +
  geom_smooth(data = combined_df_biased, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = pred_mean_b, aes(x = n_occ, y = iperv_mean), 
            color = "darkgreen", size = 1.2) +
  labs(title = "Mean (biased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()


## Plot: unbiased & biased
combined_df$total <- "unbiased"
combined_df_biased$total <- "biased"
combined_data <- rbind(combined_df, combined_df_biased)

# Filter NA
pred_mean <- pred_mean %>% filter(!is.na(n_occ) & !is.na(iperv_mean))
pred_mean_b <- pred_mean_b %>% filter(!is.na(n_occ) & !is.na(iperv_mean))


ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  
  geom_smooth(data = combined_df_biased, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  
  geom_line(data = pred_mean, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  
  geom_line(data = pred_mean_b, aes(x = n_occ, y = iperv_mean), 
            color = "darkgreen", size = 1.2) +
  
  labs(title = "Hypervolume (Unbiased vs Biased)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()

```
``` r

## Useful values 
# Convergence values
final_unbiased <- tail(pred_mean, 1)
final_biased <- tail(pred_mean_b, 1)
final_results <- data.frame(
  type = c("unbiased", "biased"),
  n_occ = c(final_unbiased$n_occ, final_biased$n_occ),
  iperv = c(final_unbiased$iperv_mean, final_biased$iperv_mean)
)


# CSV
write.csv(final_results, file = "specie_1.csv", row.names = FALSE)
```

``` r


############ Estimating the Area Of Applicability of spatial prediction models ###############
# https://hannameyer.github.io/CAST/articles/cast02-AOA-tutorial.html


### AOA for Spatially Clustered Data: Null Model vs Biased
###  Null model


## Model training
# A machine learning algorithm will be applied to learn the relationships between predictors and response

## Train data: must be converted in the format required by terra::extract
pa_points <- presence.points$sample.points[,-(3:4)] %>% as.data.frame() %>% st_as_sf(., coords = c("x","y"), crs = 4326)

# Raster data
mydata_aoa <- rast(mydata_backup)

# Subset of the original 200 points
pa_points <- pa_points[rownames(occurrences_values), ]

# From raster, extract corresponding values 
trainDat_null <- terra::extract(mydata_aoa, pa_points, na.rm = FALSE)

# From raster, extract suitability values, NA omit, assign spatial reference
trainDat_null$response <- terra::extract(random.sp$suitab.raster, pa_points, na.rm=FALSE, ID=FALSE)
trainDat_null <- data.frame(trainDat_null, pa_points) %>% na.omit()

## Train model for Spatially Clustered Data
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
folds_null <- CreateSpacetimeFolds(trainDat_null, spacevar = "geometry", k = 4)

set.seed(15)
model_null <- train(trainDat_null[,names(mydata_aoa)],
                    trainDat_null$response$`VSP suitability`,
                    method = "rf",
                    importance = TRUE,
                    tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                    trControl = trainControl(method ="cv", index = folds_null$index))

print(model_null)

# Variable Importance of each predictor
plot(varImp(model_null, scale = F), col="black", main = "Importance of each predictor", axes =FALSE)


## Predict and calculate error 
# The trained model is then used to make predictions for the entire area of interest
prediction_null <- predict(mydata_aoa, model_null, na.rm=T)

# The difference bewteen prediction and reference is the true prediction error 
truediff_null <- abs(prediction_null - random.sp$suitab.raster)

# Plot Prediction, Reference and Difference
par(mfrow = c(1, 3)) 
plot(prediction_null, main = "Prediction with RF", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(random.sp$suitab.raster, main = "Reference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))
plot(truediff_null, main = "Difference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))


dev.off()

## The AOA calculation takes the model as input to extract the importance of the predictors, 
# used as weights in multidimensional distance calculation.
AOA_null <- aoa(mydata_aoa, model_null, LPD = TRUE, verbose = FALSE)

# Features: DI, LPD, AOA
print(AOA_null)


par(mfrow = c(1, 3)) 
# DI: normalized and weighted minimum distance to a nearest training data point 
# divided by the average distance within the training data
plot(AOA_null$DI, col = viridis(100), main = "DI")

# LPD: absolute count of training data points
plot(AOA_null$LPD, col = viridis(100), main = "LPD")

# AOA: derived from the DI by using a threshold.
plot(prediction_null, col=inferno(100), main = "Prediction for Area of Applicability")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))

dev.off()


## Biased points
###############################################################################################
biased_sp_points <- points_biased %>% st_as_sf(., crs = 4326)
biased_sp_points <- biased_sp_points[,-(1:8)]

# From raster, extract corresponding values 
trainDat_biased <- terra::extract(mydata_aoa, biased_sp_points, na.rm=FALSE)

# From raster, extract suitability values 
trainDat_biased$response <- terra::extract(random.sp$suitab.raster, biased_sp_points, na.rm = FALSE, ID=FALSE)
trainDat_biased <- data.frame(trainDat_biased, biased_sp_points)

# Omit NA
trainDat_biased <- na.omit(trainDat_biased)
trainDat_biased
## Train model
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
folds_biased <- CreateSpacetimeFolds(trainDat_biased, spacevar = "geometry", k = 10)
set.seed(15)
model_biased <- train(trainDat_biased[,names(mydata_aoa)],
                      trainDat_biased$response$`VSP suitability`,
                      method = "rf",
                      importance = TRUE,
                      tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                      trControl = trainControl(method ="cv", index = folds_biased$index))

print(model_null)

# Variable Importance of each predictor
plot(varImp(model_biased, scale = F), col="black")
plotResponse(random.sp)

## Predict and calculate error 
# The trained model is then used to make predictions for the entire area of interest
prediction_biased <- predict(mydata_aoa, model_biased, na.rm=T)

# Difference bewteen prediction and reference: true prediction error 
truediff_biased <- abs(prediction_biased - random.sp$suitab.raster)

# Plot Prediction, Reference and Difference
par(mfrow = c(1, 2)) 
plot(prediction_biased, main = "Prediction with RF", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(random.sp$suitab.raster, main = "Reference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

dev.off()


# Plot Prediction, Reference and Difference
par(mfrow = c(1, 2)) 
plot(prediction_biased, main = "Prediction with RF", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(random.sp$suitab.raster, main = "Reference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))
plot(truediff_biased, main = "Difference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

par(mfrow = c(1, 2)) 
plot(prediction_null, main = "RF Null Model", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(prediction_biased, main = "RF Biased data", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))
## The AOA calculation takes the model as input to extract the importance of the predictors 
# used as weights in multidimensional distance calculation.
AOA_biased <- aoa(mydata_aoa, model_biased, LPD = TRUE, verbose = FALSE)

# Features: DI, LPD, AOA
class(AOA_biased)
print(AOA_biased)

# Plotting the aoa object 
# Shows the distribution of DI values within the training data and the DI of the new data.
plot(AOA_biased)

dev.off()

plot(truediff_biased, col = viridis(100), main = "True Prediction Error")

# DI: normalized and weighted minimum distance to a nearest training data point 
# divided by the average distance within the training data
plot(AOA_biased$DI, col = viridis(100), main = "DI")

# LPD: absolute count of training data points
plot(AOA_biased$LPD, col = viridis(100), main = "LPD")


# AOA: derived from the DI by using a threshold.
plot(prediction_biased, col=inferno(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))


dev.off()


############################ Comparison ###############################################
## Set same scale
par(mfrow=c(1,2))
plot(prediction_null, col=viridis(100), main = "Prediction for AOA (Null)")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))

plot(prediction_biased, col=viridis(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))

#######################################################################################
model_null$results
model_biased$results

## Difference? Show in the map (calc. pixel)
plot(prediction_biased)
plot(AOA_biased$AOA)

# Uniased Masked 
masked_raster_null <- mask(prediction_null, AOA_null$AOA, maskvalues=0, updatevalue=NA)

# Biased Masked
masked_raster_biased <- mask(prediction_biased, AOA_biased$AOA, maskvalues=0, updatevalue=NA)

par(mfrow=c(1,3))
plot(masked_raster_null)
plot(masked_raster_biased)

dev.off()

# Pixels in Null Model only
diff_null_only <- ifel(!is.na(masked_raster_null) & is.na(masked_raster_biased), 1, NA)

# Pixels in Biased Model only
diff_biased_only <- ifel(is.na(masked_raster_null) & !is.na(masked_raster_biased), -1, NA)

# Merge
diff_raster <- merge(diff_null_only, diff_biased_only)

# Palette
col_palette <- c("deeppink", "darkgreen")

# Plot
par(mfrow = c(1, 3), mar = c(5, 4, 4, 4) + 0.1)
plot(masked_raster_null, main = "Null", col=viridis(100), legend =FALSE)
plot(masked_raster_biased, main = "Biased", col=viridis(100), legend = FALSE)
plot(diff_raster, col = col_palette, main = "Difference", legend = FALSE)
par(mar = c(5, 4, 4, 4) + 0.1, xpd = TRUE)
legend("topleft", legend = c("Bias - Null", "Null - Bias"), fill = col_palette, cex = 0.8, bty = "n")
par(mfrow = c(1, 1))
dev.off()


#### Spatial difference ####
pixel_values <- values(diff_raster)

# Num. red and blue pixels
num_red_pixels <- sum(pixel_values == -1, na.rm = TRUE)
num_blue_pixels <- sum(pixel_values == 1, na.rm = TRUE)

# Print
cat("N. red pixels:", num_red_pixels)
cat("Area red pixels (km^2):", area_red_km2)
cat("N. blu pixels", num_blue_pixels)
cat("Area blu pixels (km^2):", area_blue_km2)
```
# References
 * Gilardi A, Lovelace R (2024). osmextract: Download and Import Open Street Map Data Extracts. R package version 0.5.1.900, [https://github.com/ropensci/osmextract](https://docs.ropensci.org/osmextract/).
