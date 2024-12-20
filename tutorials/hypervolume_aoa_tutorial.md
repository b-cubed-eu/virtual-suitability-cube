# Evaluating the impact of sampling bias on the measurement of the ecological niche
To enhance our understanding of biodiversity under changing climatic conditions, **Species Distribution Models (SDMs)**, also known as Ecological Niche Models (ENMs), have emerged as a powerful tool to characterize species' niches and quantify their ecological requirements.  

In spite of the rapid development of methodologies, theoretical constraints still exist in SDMs. For example, a significant methodological assumption required by the use of SDMs is the need for unbiased data. 

**Sampling bias** is a considerable issue when creating large-scale distribution models. It is likely to occur when data are collected without using a planned sample strategy.
One of the most recognized forms of bias in distributional data is the high concentration of observations (or collection sites) along **roads**.

The aim of this tutorial is to propose a method that shows how subsampling due to roadside bias affects niche completeness.

This will be achieved by generating a set of virtual occurrences using the R package `virtualspecies` ([Leroy, 2016](https://borisleroy.com/files/virtualspecies-tutorial.html)), simulating an erroneous dataframe with bias with respect to a street grid, and for both dataframes, calculating the **hypervolume** as occurrences increase. 

In this way, we will have two accumulation curves that relate the number of occurrences to the hypervolume for both the unbiased dataset and the biased one. 
Quantifying the hypervolume difference provides insight into the completeness of the ecological niche, but it lacks spatial information. 

To provide this information, we will use the Area Of Applicability ([Meyer, 2021](https://hannameyer.github.io/CAST/articles/cast02-AOA-tutorial.html), defined as the region where a predictive model can reliably be applied based on the relationships it learned from the training data.

The Area Of Applicability was developed to evaluate the model's performance and its ability to generalize to new areas. However, since our focus is on the data and their biases rather than the model itself, we will generate the same model using the two datasets we created before, one unbiased and one biased.

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
Let's load the vector files: the study area and the road network from Open Street Map ([osmextract R package](https://docs.ropensci.org/osmextract/)).

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
Generating random species from known environmental data allows controlling the factors that can influence the distribution of real data. This can be done through `virtualspecies` R package.  
To create a series of occurrence points for a species, it is necessary to go through 3 steps. Within each step, you can adjust parameters whose meaning is explained in the virtualspecies tutorial.

1. By intersecting bioclimatic data, the first output obtained is a **suitability map**
2. The suitability map is converted into a binary **presence/absence map** through a probability function (logistic curve) that associates the suitability value with the probability of finding the virtual species for each pixel. This subset of the environmental niche that is actually occupied by the species corresponds to the realized niche (Hutchinson, 1957).
3. The third step consists of generating, within the presence/absence raster, a series of **occurrence points**. 


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


## Ecological niche as hypervolume
Hutchinson defined an ecological niche as an n-dimensional volume in the environmental space where a species can maintain a viable population and persist along
time. The Hutchinsonian niche, or n-dimensional environmental space, is defined as the hypervolume (Blonder et al., 2014). 

The calculation of the hypervolume is performed using the `hypervolume_gaussian` function from the R package `hypervolume` [Blonder, 2014](https://onlinelibrary.wiley.com/doi/10.1111/geb.12146). `hypervolume_gaussian` constructs a hypervolume by building a gaussian kernel density estimate on an adaptive grid of random points wrapped around the original data points. 

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
# useless columns 
drops <- c("Real","Observed", "x", "y")
occurrences_values <- filtered_occ[ , !(names(filtered_occ) %in% drops)]

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
  <img width="430" height="300" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/distance_poster.png">
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
  <img width="430" height="300" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/sampling_prob_poster.png">
</p>

Ora associamo a ciascun punto la propria probabilità di essere campionato e filtriamo quelli la cui probabilità è 100%: dal campionamento randomico, si ottiene un sottoinsieme che presenta un bias spaziale, dato che è stato campionato preferenzialmente lungo il reticolo stradale

``` r
## occurrences as points
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")


# add probability value
points_biased <- coord_occ %>%
  cbind(terra::extract(prob_raster, ., ID = FALSE)) %>%
  subset(.$layer == 1)


## hypervolume of occurrences (random sampled: null model) 
# num. simulations each species
num_sim <- 3

# set stop point according to the number of biased occurrences: same sampling effort
nrow(points_biased)
nrow(occurrences_values)
stop <-  ceiling(nrow(points_biased) + 0.2 * (nrow(points_biased)))

# random subsample of occurrences from null model: 20%
occurrences_values <- occurrences_values[sample(nrow(occurrences_values), stop), ]

## plot: map with unbiased-biased points
# index
indices <- rownames(occurrences_values)
indices <- as.numeric(indices)
filtered_coord_occ <- coord_occ[indices, ]

par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 

plot(prob_r, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
points(filtered_coord_occ, cex = 0.6)
points(points_biased, col = "red", cex = 0.6)
legend("topright", legend = c("Unbiased", "Biased"), col = c("black", "red"), pch = 19, cex = 0.8,
       xpd = TRUE, y.intersp = 0.8)
```

<p align="center">
  <img width="430" height="300" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/unbiased_biased_poster.png">
</p>


## Unbiased vs biased hypervolume
We construct hypervolumes starting from 40 occurrences, increasing by 30 occurrences each time up to the maximum number.

Using the Loess method, a powerful but simple strategy for fitting smooth curves to empirical data within the `geom_smooth()` function of ggplot, we obtain as output an accumulation curve. 

If the curve saturates around a certain value, the occurrences adequately represent all suitable environments for the species (Arlè). 

Alternatively, if the species shows a more linear relationship, gaps in the spatial data might have prevented a complete representation of the observed native niche breadth (Hortal et al., 2008).  

The accumulation curve will be simulated 10 times in order to obtain an average. 

``` r
# list with the occurrences we want to test
hyp_steps <- c(seq(from = 40, to = stop, by = 30), stop)

# empty list 
all_sim <- list()

# for cycle for simulations
for (sim in 1:num_sim) {
  
  list_output <- list()
  
  for (i in seq_along(hyp_steps)) {
    d_hyp <- acc_curve(occurrences_values, hyp_steps[i])
    list_output[[i]] <- d_hyp[[1]]
    }
  
  all_sim[[sim]] <- list_output
  
}

# all simulations in one df
combined_df <- do.call(rbind, lapply(seq_along(all_sim), function(sim) {
  do.call(rbind, lapply(all_sim[[sim]], function(df) {
    df$sim <- sim
    df
  }))
}))


# mean predictions (LOESS): x sequence 
x_seq <- seq(min(combined_df$n_occ), max(combined_df$n_occ), length.out = 100)

# mean predictions: LOESS method
loess_predictions <- lapply(unique(combined_df$n_occ), function(n) {
  preds <- sapply(all_sim, function(lista) {
    loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
    predict(loess_fit, newdata = data.frame(n_occ = n))
  })
  
  data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
  
})


# mean df
pred_mean <- do.call(rbind, loess_predictions)

# plot: unbiased hypervolume
ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = pred_mean, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  labs(title = "Mean (unbiased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()

```

<p align="center">
  <img width="450" height="350" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/unbiased_poster.png">
</p>

Ora calcoliamo l'ipervolume del sottocampione biased.
``` r
## hypervolume of biased occurrences (road driven: biased sampling)
biased_df <- points_biased %>%
  as.data.frame() %>%
  .[,-c(5:8)]

# stop
stop_biased <- nrow(biased_df)
hyp_steps_b <- c(seq(from = 20, to = stop_biased, by = 20), stop_biased)

# empty list
all_sim_b <- list()

# for cycle for simulations
for (sim in 1:num_sim) {
  
  list_output_b <- list()
  
  for (i in seq_along(hyp_steps_b)) {
    d_hyp <- acc_curve(biased_df, hyp_steps_b[i])
    list_output_b[[i]] <- d_hyp[[1]]
    }
  
  all_sim_b[[sim]] <- list_output_b
  
}

# combined df
combined_df_biased <- do.call(rbind, lapply(seq_along(all_sim_b), function(sim) {
  
  do.call(rbind, lapply(all_sim_b[[sim]], function(df) {
    df$sim <- sim
    df
    
  }))
  
}))


# mean LOESS
loess_predictions_biased <- lapply(unique(combined_df_biased$n_occ), function(n) {
  
  preds <- sapply(all_sim_b, function(lista) {
    loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
    predict(loess_fit, newdata = data.frame(n_occ = n))
  })
  
  data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
  
})

# mean in one df
pred_mean_b <- do.call(rbind, loess_predictions_biased)

## plot: biased hypervolume
ggplot() +
  geom_smooth(data = combined_df_biased, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = pred_mean_b, aes(x = n_occ, y = iperv_mean), 
            color = "darkgreen", size = 1.2) +
  labs(title = "Mean (biased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()
```
<p align="center">
  <img width="450" height="350" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/biased_poster.png">
</p>

Mettiamo nello stesso grafico i due ipervolumi, ricordando che fanno riferimento alla stessa specie
``` r
# plot: unbiased & biased
combined_df$total <- "unbiased"
combined_df_biased$total <- "biased"
combined_data <- rbind(combined_df, combined_df_biased)

# filter NA
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

<p align="center">
  <img width="500" height="350" src="https://github.com/rociobeatrizc/virtual-suitability-cube/blob/main/images/2_hyperv_poster.png">
</p>

# Area of Applicability
We want to utilize this method not to evaluate the model's performance but to assess the quality of the starting data based on their predictive capacity when used to build a model.

We used Random Forest (Breiman, 2001) as our machine learning algorithm because it is widely used in environmental mapping (e.g., Bastin et al., 2019). For model training, we extracted the 4 predictors from the sampling data points.

We obtain a subset of the original study area where suitability estimates can be made. By repeating this process for both the unbiased and the biased datasets, we can calculate the difference in the number of pixels between the two outputs. We expect that the area of applicability of a model calibrated on a randomly sampled dataset will be larger compared to that of a biased dataset

The model built on the unbiased dataset will be called the **null model**, while the one built on the biased dataset will be the **biased model**.

``` r
### AOA for Spatially Clustered Data: Null Model vs Biased
###  null model

## model training
# a machine learning algorithm will be applied to learn the relationships between predictors and response

## train data: must be converted in the format required by terra::extract
pa_points <- presence.points$sample.points[,-(3:4)] %>% as.data.frame() %>% st_as_sf(., coords = c("x","y"), crs = 4326)

# raster data
mydata_aoa <- rast(mydata_backup)

# subset of the original 200 points
pa_points <- pa_points[rownames(occurrences_values), ]

# from raster, extract corresponding values 
trainDat_null <- terra::extract(mydata_aoa, pa_points, na.rm = FALSE)

# from raster, extract suitability values, NA omit, assign spatial reference
trainDat_null$response <- terra::extract(random.sp$suitab.raster, pa_points, na.rm=FALSE, ID=FALSE)
trainDat_null <- data.frame(trainDat_null, pa_points) %>% na.omit()

## train model for Spatially Clustered Data
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
folds_null <- CreateSpacetimeFolds(trainDat_null, spacevar = "geometry", k = 4)

set.seed(15)
model_null <- train(trainDat_null[,names(mydata_aoa)],
                    trainDat_null$response$`VSP suitability`,
                    method = "rf",
                    importance = TRUE,
                    tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                    trControl = trainControl(method ="cv", index = folds_null$index))

## predict and calculate error 
# the trained model is then used to make predictions for the entire area of interest
prediction_null <- predict(mydata_aoa, model_null, na.rm=T)

# the difference bewteen prediction and reference is the true prediction error 
truediff_null <- abs(prediction_null - random.sp$suitab.raster)

## the AOA calculation takes the model as input to extract the importance of the predictors 
# used as weights in multidimensional distance calculation.
AOA_null <- aoa(mydata_aoa, model_null, LPD = TRUE, verbose = FALSE)

# AOA: derived from the DI by using a threshold.
plot(prediction_null, col=inferno(100), main = "Prediction for Area of Applicability")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))
```

<p align="center">
  <img width="500" height="450" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/aoa_unbiased_poster.png">
</p>

Thanks to the AOA, we can see the extent of the area where Random Forest is capable of making predictions
Now, let's see how it works if we train the algorithm on biased data

``` r
## biased points
biased_sp_points <- points_biased %>% st_as_sf(., crs = 4326)
biased_sp_points <- biased_sp_points[,-(1:8)]

# from raster, extract corresponding values 
trainDat_biased <- terra::extract(mydata_aoa, biased_sp_points, na.rm=FALSE)

# from raster, extract suitability values 
trainDat_biased$response <- terra::extract(random.sp$suitab.raster, biased_sp_points, na.rm = FALSE, ID=FALSE)
trainDat_biased <- data.frame(trainDat_biased, biased_sp_points)

# omit NA
trainDat_biased <- na.omit(trainDat_biased)

## train model
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
folds_biased <- CreateSpacetimeFolds(trainDat_biased, spacevar = "geometry", k = 10)
set.seed(15)
model_biased <- train(trainDat_biased[,names(mydata_aoa)],
                      trainDat_biased$response$`VSP suitability`,
                      method = "rf",
                      importance = TRUE,
                      tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                      trControl = trainControl(method ="cv", index = folds_biased$index))

## predict and calculate error 
# the trained model is then used to make predictions for the entire area of interest
prediction_biased <- predict(mydata_aoa, model_biased, na.rm=T)

# difference bewteen prediction and reference: true prediction error 
truediff_biased <- abs(prediction_biased - random.sp$suitab.raster)

# Random Forest trained on unbiased and biased data
par(mfrow = c(1, 2)) 
plot(prediction_null, main = "RF Null Model", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(prediction_biased, main = "RF Biased data", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

```

<p align="center">
  <img width="600" height="400" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/rf_null_biased_poster.png">
</p>


``` r

## the AOA calculation takes the model as input to extract the importance of the predictors 
# used as weights in multidimensional distance calculation.
AOA_biased <- aoa(mydata_aoa, model_biased, LPD = TRUE, verbose = FALSE)

# AOA: derived from the DI by using a threshold.
plot(prediction_biased, col=inferno(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))
```

<p align="center">
  <img width="500" height="350" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/aoa_biased_poster.png">
</p>

If we compare the two outputs, we can immediately see that there is a difference in the extention of the prediction area: 

```  r

par(mfrow=c(1,2))
plot(prediction_null, col=viridis(100), main = "Prediction for AOA (Null)")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))

plot(prediction_biased, col=viridis(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))
```
<p align="center">
  <img width="600" height="400" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/aoa_biased_unbiased_poster.png">
</p>

Finally, let's see the different amount of pixels 
``` r

## difference? Show in the map (calc. pixel)
# unbiased masked 
masked_raster_null <- mask(prediction_null, AOA_null$AOA, maskvalues=0, updatevalue=NA)

# biased masked
masked_raster_biased <- mask(prediction_biased, AOA_biased$AOA, maskvalues=0, updatevalue=NA)

# pixels in null model only
diff_null_only <- ifel(!is.na(masked_raster_null) & is.na(masked_raster_biased), 1, NA)

# pixels in biased model only
diff_biased_only <- ifel(is.na(masked_raster_null) & !is.na(masked_raster_biased), -1, NA)

# merge
diff_raster <- merge(diff_null_only, diff_biased_only)

# palette
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
```

<p align="center">
  <img width="600" height="400" src="https://github.com/b-cubed-eu/virtual-suitability-cube/blob/main/images/difference.png">
</p>
