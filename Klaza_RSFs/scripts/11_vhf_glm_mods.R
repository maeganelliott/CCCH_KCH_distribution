#VHF GLM model selection - Klaza herd
#Oct 24, 2024
#Reviewed Nov 20, 2024

library(easystats)
library(see)
library(effectsize)
library(parameters)
library(plyr)
library(tidyverse)
library(terra)
library(sf)
library(dplyr)
library(sp)
library(corrplot)
library(report)
library(terra)
library(predicts)
library(ggplot2)
library(car)
library(AICcmodavg)

set.seed(123)

#############################################################################
##SNOW SEASON

#Import used and available points with covariates previously extracted to them
snow_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_snow_data.csv')
#check data:
str(snow_data) #numeric and integers, will work with scale function
snow_data$Band_1 <- as.numeric(snow_data$Band_1) #set elevation as numeric
colSums(is.na(snow_data)) #check for NAs - 2 points with NAs to remove. 
snow_data <- na.omit(snow_data)


#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(snow_data), 0.8*nrow(snow_data))
snow_train = snow_data[n,]
snow_test = snow_data[-n,] 

colSums(is.na(snow_train))
colSums(is.na(snow_test))

##########################################################################
#Run univariate models to determine best scale for each PFT veg raster

#Lichen
lichbi30 <- glm(presence~lichen_bi95, family=binomial, data=snow_train)
lichenbi100 <- glm(presence~lichen_bi95_100, family=binomial, data=snow_train)
lichenbi250 <- glm(presence~lichen_bi95_250, family=binomial, data=snow_train)
summary(lichenbi250) # Lichen bivariate 250 m is the best based on AIC model selection. Medium, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi95, family=binomial, data=snow_train)
conbi100 <- glm(presence~conifer_bi95_100, family=binomial, data=snow_train)
conbi250 <- glm(presence~conifer_bi95_250, family=binomial, data=snow_train)
summary(conbi30) # Conifer bivariate 30 m is the best based on AIC model selection but all scales very similar. Positive effect.  
#AIC 
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~decid_bi95, family=binomial, data=snow_train)
decidbi100 <- glm(presence~decid_bi95_100, family=binomial, data=snow_train)
decidbi250 <- glm(presence~decid_bi95_250, family=binomial, data=snow_train)
summary(decidbi30) # decid bivariate 30 m is the best based on AIC model selection. Positive effect.  Exclude from snow model, not biologically relevant.. 
#AIC 
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#Graminoid
grambi30 <- glm(presence~gram_bi95, family=binomial, data=snow_train)
grambi100 <- glm(presence~gram_bi95_100, family=binomial, data=snow_train)
grambi250 <- glm(presence~gram_bi95_250, family=binomial, data=snow_train)
summary(grambi250) # Graminoid bivariate 250 m is the best based on AIC model selection but all scales very similar. Positive effect.  
#AIC 
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glms_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glms_names)

#Include: dem, slope, lichenbi250, coniferbi30, grambi250, fire
#exclude deciduous shrub for winter?

#Check Correlation of Snow season variables
#add quadratic terms
snow_data$slope2 <- snow_data$slope^2 #quadratic transformation
snow_data$dem2 <- snow_data$Band_1^2 #remember elevation is called Band_1

missing_cols <- setdiff(c("Band_1", "dem2", "slope", "slope2", "lichen_bi95_250", "conifer_bi95", "decid_bi95", "gram_bi95_250", "fire1995", "dist1995_30", "dist1995_500", "dist1995_1000", "dist1995_2000", "dist1995_3000", "lineden1995_250", "lineden1995_500", "lineden1995_1000", "lineden1995_2000", "lineden1995_3000"), colnames(snow_data))
missing_cols # check for miss-matched names if correlation code doesn't work.


snow_dat.cor = snow_data[,c("Band_1", "dem2", "slope", "slope2", "lichen_bi95_250", "conifer_bi95", "decid_bi95", "gram_bi95_250", "fire1995", "dist1995_30", "dist1995_500", "dist1995_1000", "dist1995_2000", "dist1995_3000", "lineden1995_250", "lineden1995_500", "lineden1995_1000", "lineden1995_2000", "lineden1995_3000")]
names(snow_dat.cor) = c("Band_1", "dem2", "slope", "slope2", "lichen_bi95_250", "conifer_bi95", "decid_bi95", "gram_bi95_250", "fire1995", "dist1995_30", "dist1995_500", "dist1995_1000", "dist1995_2000", "dist1995_3000", "lineden1995_250", "lineden1995_500", "lineden1995_1000", "lineden1995_2000", "lineden1995_3000")
cor_matrix <- cor(snow_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(snow_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/stats/vhf_snow_correlation_matrix.csv", row.names = TRUE)
#some correlation between disturbance variables to watch for. 

############################################################################
#Snow Season GLM Models

#Some experimentation. 
# all variables - dem, slope, lichenbi250, coniferbi30, grambi250, fire
s_test1 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(I(slope^2)) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#remove elevation
s_test2 <- glm(presence~scale(slope) + scale(I(slope^2)) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#remove quadratic slope term
s_test3 <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#lichen on its own
s_test4 <- glm(presence~scale(lichen_bi95_250), family=binomial,data=snow_train)
#elevation on its own
s_test5 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)), family=binomial,data=snow_train)

summary(s_test1)

#In full model lichen has a positive effect and is close to significant (0.05), conifer has a significant positive effect, graminoid has 
#significant positive effect, fire has significant negative effect. Linear slope term is negative and close to significant (0.09). 
#When quadratic term is removed linear slope term is significant. Best depicted by a linear relationship. 
#Conifer cover and slope seem to be main drivers of model, followed by graminoid and fire. 
#lichen modeled on its own has a significant positive effect
#elevation modeled on its own has a significant quadratic relationship. 


#GLM - based on covariate exploration / distributions
#Include: elevation (q), slope, lichen, conifer, graminoid, fire.
#evaluate with AIC model selection - see if elevation improves model fit

#global 
snow_glm_1e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#primary model drivers (slope, conifer, graminoid, fire); no elevation and lichen.
snow_glm_2e <- glm(presence~scale(slope) + scale(conifer_bi95) + scale(gram_bi95_100) + scale(fire1995), family=binomial,data=snow_train)
#all variables minus fire
snow_glm_3e  <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250), family=binomial,data=snow_train)
#all variables minus graminoid
snow_glm_4e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(fire1995), family=binomial,data=snow_train)
#all variables minus conifer
snow_glm_5e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_250) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#all variables minus lichen
snow_glm_6e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#all variables minus slope
snow_glm_7e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#all variables minus elevation
snow_glm_8e <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#null
snow_glm_9e <- glm(presence~1, family=binomial,data=snow_train)

#AIC Model Selection
snow_e_glms <- list(snow_glm_1e, snow_glm_2e, snow_glm_3e, snow_glm_4e, snow_glm_5e, snow_glm_6e, snow_glm_7e, snow_glm_8e, snow_glm_9e)
snow_e_glms_names <- c('snow_glm_1e', 'snow_glm_2e', 'snow_glm_3e', 'snow_glm_4e', 'snow_glm_5e', 'snow_glm_6e', 'snow_glm_7e', 'snow_glm_8e', 'snow_glm_9e')
aictab(cand.set = snow_e_glms, modnames = snow_e_glms_names)

summary(snow_glm_8e)
#removal of elevation improves model fit (8e)
#global model is second best
#model without lichen is third best
#fire, slope, graminoid and conifer are most important variables.

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(snow_glm_8e, newdata = snow_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = snow_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

# Step 4: Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(snow_glm_8e, newdata = snow_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = snow_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Step 4: Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()

#model tends to underpredict due to small sample size. Conservative. 

#Evaluation with AUC / ROC
#unscaled model:
snow_glm_8e_us <- glm(presence~slope + lichen_bi95_250 + conifer_bi95 + gram_bi95_250 + fire1995, family=binomial,data=snow_train)
#evaluate
snow_glm_8e_eval <- pa_evaluate(predict(snow_glm_8e_us, snow_test[snow_test$presence==1, ]), predict(snow_glm_8e_us, snow_test[snow_test$presence==0, ]))
print(snow_glm_8e_eval)
plot(snow_glm_8e_eval, "ROC")
#AUC=0.701 #moderate discrimination / fair model performance. 
#Model performance with easystats
r2(snow_glm_8e_us) #Tjur's R2 = 0.067

#Caribou occurrence prediction map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/k_vhf_mask_30m_2.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters1995_30m.tif')
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters1995_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
slope <- rast30m$slope
lichen_bi95_250 <- rast250m$lichen_bi95_250 |>
  resample(mask, method='near')
conifer_bi95 <- rast30m$conifer_bi95
gram_bi95_250 <- rast250m$gram_bi95_250 |>
  resample(mask, method='near')
fire1995 <- rast30m$fire1995

rasters <- c(slope, lichen_bi95_250, conifer_bi95, gram_bi95_250, fire1995) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, snow_glm_8e_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_snow_mod8e_prediction.tif')

#Add gps points to map
#filtered used points w/o covariates extracted to them:
snow_points <- read.csv('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/k_vhf_snow.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snow_points)
# convert from lat/long to UTM to match rasters
snow_points <- st_as_sf(snow_points, coords = c("location.long", "location.lat"), crs = 4326)
snow_points <- st_transform(snow_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snow_points)
plot(st_geometry(snow_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

###Phase 2 - Disturbance Model Selection####

#test effect of disturbance variables on base GLM 8e. 
#In combined models ensure disturbance variables are not correlated. 

#glm 8e + buffered dist 30 m
snow_glm_1d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(dist1995_30), family=binomial,data=snow_train)
#glm 8e + buffered dist 500 m
snow_glm_2d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(dist1995_500), family=binomial,data=snow_train)
#glm 8e + buffered dist 1000 m
snow_glm_3d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(dist1995_1000), family=binomial,data=snow_train)
#glm 8e + buffered dist 2000 m
snow_glm_4d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(dist1995_2000), family=binomial,data=snow_train)
#glm 8e + buffered dist 3000 m
snow_glm_5d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(dist1995_3000), family=binomial,data=snow_train)
#glm 8e + lineden 250 m
snow_glm_6d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(lineden1995_250), family=binomial,data=snow_train)
#glm 8e + lineden 500 m
snow_glm_7d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(lineden1995_500), family=binomial,data=snow_train)
#glm 8e + lineden 1000 m 
snow_glm_8d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(lineden1995_1000), family=binomial,data=snow_train)
#glm 8e + lineden 2000 m
snow_glm_9d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(lineden1995_2000), family=binomial,data=snow_train)
#glm 8e + lineden 3000 m
snow_glm_10d <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995) + scale(lineden1995_3000), family=binomial,data=snow_train)
#glm 8e
snow_glm_8e <- glm(presence~scale(slope) + scale(lichen_bi95_250) + scale(conifer_bi95) + scale(gram_bi95_250) + scale(fire1995), family=binomial,data=snow_train)
#null
snow_glm_11d <- glm(presence~1, family=binomial,data=snow_train)

snow_d_glms <- list(snow_glm_1d, snow_glm_2d, snow_glm_3d, snow_glm_4d, snow_glm_5d, snow_glm_6d, snow_glm_7d, snow_glm_8d, snow_glm_9d, snow_glm_10d, snow_glm_11d, snow_glm_8e)
snow_d_glms_names <- c('snow_glm_1d', 'snow_glm_2d', 'snow_glm_3d', 'snow_glm_4d', 'snow_glm_5d', 'snow_glm_6d', 'snow_glm_7d', 'snow_glm_8d', 'snow_glm_9d', 'snow_glm_10d', 'snow_glm_11d', 'snow_glm_8e')
aictab(cand.set = snow_d_glms, modnames = snow_d_glms_names)

summary(snow_glm_3d)
#lineden 2000 m radii in top model followed by lineden 3000 m radii, both negative but non-significant effect. 
#buffered disturbance 1000 m in third best model, significant negative effect
#buffered disturbance 30 m in fourth best model, negative but not significant effect. 
#buffered disturbance 1000 m is correlated with lineden 2000 and 3000 m so can't go in combined model. 
#proceed with GLM 3d

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(snow_glm_3d, newdata = snow_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = snow_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

# Step 4: Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(snow_glm_3d, newdata = snow_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = snow_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Step 4: Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()

#model tends to underpredict due to small sample size. Conservative.

#Evaluation with AUC / ROC
#unscaled model:
snow_glm_3d_us <- glm(presence~slope + lichen_bi95_250 + conifer_bi95 + gram_bi95_250 + fire1995 + dist1995_1000, family=binomial,data=snow_train)
#evaluate
snow_glm_3d_eval <- pa_evaluate(predict(snow_glm_3d_us, snow_test[snow_test$presence==1, ]), predict(snow_glm_3d_us, snow_test[snow_test$presence==0, ]))
print(snow_glm_3d_eval)
plot(snow_glm_3d_eval, "ROC")
#AUC=0.705 #moderate discrimination / fair model performance. Improved from base model. 
#Model performance and plots with easystats
r2(snow_glm_3d_us) #Tjur's R2 = 0.071

windows()
check_model(snow_glm_3d_us)
#effect size plot with see pkg from easystats
plot(effectsize(snow_glm_3d_us)) +
  ggplot2::labs(title = "Klaza - Snow 1995") +
  ggplot2::scale_y_discrete(labels = c(
    "slope" = "Slope",
    "lichen_bi95_250" = "Lichen 250 m",
    "conifer_bi95" = "Conifer 30 m",
    "gram_bi95_250" = "Graminoid 250 m",
    "fire1995" = "Burns ≤ 50 yrs",
    "dist1995_1000" = "Disturbance 1000 m"
  )) 

plot(parameters(snow_glm_3d_us)) +
  ggplot2::labs(title = "Snow")

#Save models
saveRDS(snow_glm_3d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_snow_scaled.rds")
saveRDS(snow_glm_3d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_snow_unscaled.rds")

#Caribou occurrence prediction map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/k_vhf_mask_30m_2.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters1995_30m.tif')
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters1995_250m2.tif') %>%
  resample(mask, method="near")
distrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kdistbuffrast1995_1.tif')

#prepare rasters
slope <- rast30m$slope
lichen_bi95_250 <- rast250m$lichen_bi95_250 |>
  resample(mask, method='near')
conifer_bi95 <- rast30m$conifer_bi95
gram_bi95_250 <- rast250m$gram_bi95_250 |>
  resample(mask, method='near')
fire1995 <- rast30m$fire1995
dist1995_1000 <- distrast$dist1995_1000

rasters <- c(slope, lichen_bi95_250, conifer_bi95, gram_bi95_250, fire1995, dist1995_1000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, snow_glm_3d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_snow_mod3d_prediction.tif')

#Add gps points to map
#filtered used points w/o covariates extracted to them:
snow_points <- read.csv('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/k_vhf_snow.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snow_points)
# convert from lat/long to UTM to match rasters
snow_points <- st_as_sf(snow_points, coords = c("location.long", "location.lat"), crs = 4326)
snow_points <- st_transform(snow_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snow_points)
plot(st_geometry(snow_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

########################################################################################
##SNOW FREE SEASON

#Import used and available points
#Import used and available points with covariates previously extracted to them
snowfree_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_snowfree_data.csv')
#check data:
str(snowfree_data) #numeric and integers, will work with scale function
snowfree_data$Band_1 <- as.numeric(snowfree_data$Band_1) #set elevation as numeric
colSums(is.na(snowfree_data)) #check for NAs - four with NAs
snowfree_data <- na.omit(snowfree_data) #remove NAs

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(snowfree_data), 0.8*nrow(snowfree_data))
sf_train = snowfree_data[n,]
sf_test = snowfree_data[-n,] 
###############################################################################
#Run univariate models to determine best scale for each PFT veg raster

#Lichen
lichbi30 <- glm(presence~lichen_bi95, family=binomial, data=snowfree_train)
lichenbi100 <- glm(presence~lichen_bi95_100, family=binomial, data=snowfree_train)
lichenbi250 <- glm(presence~lichen_bi95_250, family=binomial, data=snowfree_train)
summary(lichenbi100) # Lichen bivariate 100 m is the best based on AIC model selection. Medium, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi95, family=binomial, data=snowfree_train)
conbi100 <- glm(presence~conifer_bi95_100, family=binomial, data=snowfree_train)
conbi250 <- glm(presence~conifer_bi95_250, family=binomial, data=snowfree_train)
summary(conbi30) # Conifer bivariate 30 m is the best based on AIC model selection but all scales very similar. Negative effect.  
#AIC 
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~decid_bi95, family=binomial, data=snowfree_train)
decidbi100 <- glm(presence~decid_bi95_100, family=binomial, data=snowfree_train)
decidbi250 <- glm(presence~decid_bi95_250, family=binomial, data=snowfree_train)
summary(decidbi30) # decid bivariate 30 m is the best based on AIC model selection. Positive effect.  
#AIC 
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#Graminoid
grambi30 <- glm(presence~gram_bi95, family=binomial, data=snowfree_train)
grambi100 <- glm(presence~gram_bi95_100, family=binomial, data=snowfree_train)
grambi250 <- glm(presence~gram_bi95_250, family=binomial, data=snowfree_train)
summary(grambi30) # Graminoid bivariate 30 m is the best based on AIC model selection but all scales very similar. Positive effect.  
#AIC 
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glms_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glms_names)

#Include: dem, slope, lichenbi100, coniferbi30, decidbi30, grambi30 
#all variable distributions are making sense
#excluding fire for snow free season

#Check Correlation of Snowfree season variables
#add quadratic terms
snowfree_data$slope2 <- snowfree_data$slope^2 #quadratic transformation
snowfree_data$dem2 <- snowfree_data$Band_1^2 #remember elevation is called Band_1

missing_cols <- setdiff(c("Band_1", "dem2", "slope", "slope2", "lichen_bi95_100", "conifer_bi95", "decid_bi95", "gram_bi95", "fire1995", "dist1995_30", "dist1995_500", "dist1995_1000", "dist1995_2000", "dist1995_3000", "lineden1995_250", "lineden1995_500", "lineden1995_1000", "lineden1995_2000", "lineden1995_3000"), colnames(snowfree_data))
missing_cols # check for miss-matched names if correlation code doesn't work.


snowfree_dat.cor = snowfree_data[,c("Band_1", "dem2", "slope", "slope2", "lichen_bi95_100", "conifer_bi95", "decid_bi95", "gram_bi95", "fire1995", "dist1995_30", "dist1995_500", "dist1995_1000", "dist1995_2000", "dist1995_3000", "lineden1995_250", "lineden1995_500", "lineden1995_1000", "lineden1995_2000", "lineden1995_3000")]
names(snow_dat.cor) = c("Band_1", "dem2", "slope", "slope2", "lichen_bi95_100", "conifer_bi95", "decid_bi95", "gram_bi95", "fire1995", "dist1995_30", "dist1995_500", "dist1995_1000", "dist1995_2000", "dist1995_3000", "lineden1995_250", "lineden1995_500", "lineden1995_1000", "lineden1995_2000", "lineden1995_3000")
cor_matrix <- cor(snowfree_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(snow_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/stats/vhf_snowfree_correlation_matrix.csv", row.names = TRUE)
#some correlation between disturbance variables to watch for. 

##########################################################################
#GLM Model Selection

#Some experimentation
#all variables - dem, slope, lichenbi100, coniferbi30, decidbi30, grambi30
sf_test1 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(I(slope^2)) + scale(lichen_bi95_100) + scale(conifer_bi95) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#remove lichen and conifer
sf_test2 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(I(slope^2)) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#elevation only
sf_test3 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)), family=binomial,data=sf_train)
#slope only
sf_test4 <- glm(presence~scale(slope) + scale(I(slope^2)), family=binomial,data=sf_train)
#elevation, decid, gram
sf_test5 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)


summary(sf_test1)

#In full model linear elevation term (positive), decid (positive) and gram (positive) are only significant variables. 
#lichen and conifer have smallest effect size in full model. elevation is primary driver of model. 
#Graminoid close to significant (0.09) (positive) when conifer and lichen removed
#both elevation terms significant in elevation only model. Positive then negative (quadratic) relationship. 
#neither slope term is significant in slope only model, and both have a positive effect. Modelling the linear slope term (no quadratic term) it is significant. Slope seems to be best modeled by the linear term only.
#all coefficients make biological sense for the snow free season, some just aren't significant. Will proceed with AIC model selection using all identified variables. 

#Environmental Base Model Selection
#Include: dem, slope, lichenbi100, coniferbi30, decidbi30, grambi30
#Evaluate using AIC model selection

#global 
sf_glm_1e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_100) + scale(conifer_bi95) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#primary model drivers (elevation, slope, decid)
sf_glm_2e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95), family=binomial,data=sf_train)
#all variables minus graminoid
sf_glm_3e  <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_100) + scale(conifer_bi95) + scale(decid_bi95), family=binomial,data=sf_train)
#all variables minus decid
sf_glm_4e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_100) + scale(conifer_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#all variables minus conifer
sf_glm_5e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(lichen_bi95_100) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#all variables minus lichen
sf_glm_6e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(conifer_bi95) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#all variables minus slope
sf_glm_7e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(lichen_bi95_100) + scale(conifer_bi95) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#all variables minus elevation
sf_glm_8e <- glm(presence~scale(slope) + scale(lichen_bi95_100) + scale(conifer_bi95) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#no lichen and no conifer (weakest predictors)
sf_glm_9e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#null
sf_glm_10e <- glm(presence~1, family=binomial,data=sf_train)

#AIC Model Selection
sf_e_glms <- list(sf_glm_1e, sf_glm_2e, sf_glm_3e, sf_glm_4e, sf_glm_5e, sf_glm_6e, sf_glm_7e, sf_glm_8e, sf_glm_9e, sf_glm_10e)
sf_e_glms_names <- c('sf_glm_1e', 'sf_glm_2e', 'sf_glm_3e', 'sf_glm_4e', 'sf_glm_5e', 'sf_glm_6e', 'sf_glm_7e', 'sf_glm_8e', 'sf_glm_9e', 'sf_glm_10e')
aictab(cand.set = sf_e_glms, modnames = sf_e_glms_names)

summary(sf_glm_9e)

#6e is best model, excluding lichen, but delta AIC comparable to other top models
#9e is close second, excluding lichen and conifer
#third best model (2e) has elevation, slope and deciduous shrub only. 
#lichen  and conifer are weakest predictors
#delta AICc comparable among top 4 models. 
#elevation, slope, decid are most important predictors.
#proced with glm 9e

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(sf_glm_9e, newdata = sf_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = sf_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

# Step 4: Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(sf_glm_9e, newdata = sf_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = sf_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Step 4: Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()

#model tends to underpredict due to small sample size. Conservative. Deviation increases in upper bins. 

#Evaluation with AUC / ROC
#unscaled model:
sf_glm_9e_us <- glm(presence~Band_1 + I(Band_1^2) + slope + decid_bi95 + gram_bi95, family=binomial,data=sf_train)
#evaluate
sf_glm_9e_eval <- pa_evaluate(predict(sf_glm_9e_us, sf_test[sf_test$presence==1, ]), predict(sf_glm_9e_us, sf_test[sf_test$presence==0, ]))
print(sf_glm_9e_eval)
plot(sf_glm_9e_eval, "ROC")
#AUC=0.7 #moderate discrimination / fair model performance. Improved from base model.
#Model performance with easystats
r2(sf_glm_9e_us) #Tjur's R2 = 0.110


#Caribou occurrence prediction map
#load rasters 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/k_vhf_mask_30m_2.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters1995_30m.tif')
kdemrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kdem_vhf_30m.tif')

#prepare rasters
Band_1 <- as.numeric(kdemrast$Band_1)
slope <- rast30m$slope
decid_bi95 <- rast30m$decid_bi95
gram_bi95 <- rast30m$gram_bi95


rasters <- c(Band_1, slope, decid_bi95, gram_bi95) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, sf_glm_9e_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_snowfree_mod9e_prediction.tif')

#Add gps points to map
#filtered used points w/o covariates extracted to them:
snowfree_points <- read.csv('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/k_vhf_snowfree.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snowfree_points)
# convert from lat/long to UTM to match rasters
snowfree_points <- st_as_sf(snowfree_points, coords = c("location.long", "location.lat"), crs = 4326)
snowfree_points <- st_transform(snowfree_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snow_points)
plot(st_geometry(snowfree_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

##Disturbance Model Selection
#double check for collinearity before putting disturbance variables in same model

#glm 9e + dist buff 30
sf_glm_1d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(dist1995_30), family=binomial,data=sf_train)
#glm 9e + dist buff 500
sf_glm_2d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(dist1995_500), family=binomial,data=sf_train)
#glm 9e + dist buff 1000
sf_glm_3d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(dist1995_1000), family=binomial,data=sf_train)
#glm 9e + dist buff 2000
sf_glm_4d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(dist1995_2000), family=binomial,data=sf_train)
#glm 9e + dist buff 3000
sf_glm_5d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(dist1995_3000), family=binomial,data=sf_train)
#glm 9e + lineden 250
sf_glm_6d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(lineden1995_250), family=binomial,data=sf_train)
#glm 9e + lineden 500
sf_glm_7d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(lineden1995_500), family=binomial,data=sf_train)
#glm 9e + lineden 1000
sf_glm_8d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(lineden1995_1000), family=binomial,data=sf_train)
#glm 9e + lineden 2000
sf_glm_9d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(lineden1995_2000), family=binomial,data=sf_train)
#glm 9e + lineden 3000
sf_glm_10d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(lineden1995_3000), family=binomial,data=sf_train)
#glm 9e + lineden 1000 + dist buff 30
sf_glm_11d <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95) + scale(lineden1995_1000) + scale(dist1995_30), family=binomial,data=sf_train)
#glm 9e
sf_glm9e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(slope) + scale(decid_bi95) + scale(gram_bi95), family=binomial,data=sf_train)
#null
sf_glm_12d <- glm(presence~1, family=binomial,data=sf_train)

#AIC Model Selection
sf_d_glms <- list(sf_glm_1d, sf_glm_2d, sf_glm_3d, sf_glm_4d, sf_glm_5d, sf_glm_6d, sf_glm_7d, sf_glm_8d, sf_glm_9d, sf_glm_10d, sf_glm_11d, sf_glm_12d, sf_glm_9e)
sf_d_glms_names <- c('sf_glm_1d', 'sf_glm_2d', 'sf_glm_3d', 'sf_glm_4d', 'sf_glm_5d', 'sf_glm_6d', 'sf_glm_7d', 'sf_glm_8d', 'sf_glm_9d', 'sf_glm_10d', 'sf_glm_11d', 'sf_glm_12d', 'sf_glm_9e')
aictab(cand.set = sf_d_glms, modnames = sf_d_glms_names)

summary(sf_glm_8d)
#line density 1000 m radii top disturbance variable, but ranked close to lineden 2000, lineden 3000 and dist buff 30
#All top performing disturbance variables have negative effects but none are significant
#most disturbance variables decreased model fit from the base environmental model
#proceed with glm 8d with lineden 1000

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(sf_glm_8d, newdata = sf_test, type = "response")

# Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = sf_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

# Step 4: Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(sf_glm_8d, newdata = sf_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = sf_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Step 4: Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()

#model tends to under predict due to small sample size. Conservative. Deviation increases in upper bins. 

#Evaluation with AUC / ROC
#unscaled model:
sf_glm_8d_us <- glm(presence~Band_1 + I(Band_1^2) + slope + decid_bi95 + gram_bi95 + lineden1995_1000, family=binomial,data=sf_train)
#evaluate
sf_glm_8d_eval <- pa_evaluate(predict(sf_glm_8d_us, sf_test[sf_test$presence==1, ]), predict(sf_glm_8d_us, sf_test[sf_test$presence==0, ]))
print(sf_glm_8d_eval)
plot(sf_glm_8d_eval, "ROC")
#AUC=0.698 #moderate discrimination. Decreased slightly from base model.
###############################
#Model performance and plots with easystats
r2(sf_glm_8d_us) #Tjur's R2 = 0.136

windows()
check_model(sf_glm_8d_us)
#effect size plot with see pkg from easystats
plot(effectsize(sf_glm_8d_us)) +
  ggplot2::labs(title = "Klaza - Snow-free 1995") +
  ggplot2::scale_y_discrete(labels = c(
    "Band_1" = "Elevation",
    "I(Band_1^2)" = "Elevation^2",
    "slope" = "Slope",
    "decid_bi95" = "Deciduous Shrub 30 m",
    "gram_bi95" = "Graminoid 30 m",
    "lineden1995_1000" = "Line Density 1000 m"
  )) 

plot(parameters(sf_glm_8d_us)) +
  ggplot2::labs(title = "Snow-free")

#Save models
saveRDS(sf_glm_8d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_snowfree_scaled.rds")
saveRDS(sf_glm_8d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_snowfree_unscaled.rds")

#Caribou occurrence prediction map
#load rasters 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/k_vhf_mask_30m_2.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters1995_30m.tif')
kdemrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kdem_vhf_30m.tif')
linedenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/klinedenrast1995_1.tif')

#prepare rasters
Band_1 <- as.numeric(kdemrast$Band_1)
slope <- rast30m$slope
decid_bi95 <- rast30m$decid_bi95
gram_bi95 <- rast30m$gram_bi95
lineden1995_1000 <- linedenrast$lineden1995_1000


rasters <- c(Band_1, slope, decid_bi95, gram_bi95, lineden1995_1000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, sf_glm_8d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_snowfree_mod8d_prediction.tif')

#Add gps points to map
#filtered used points w/o covariates extracted to them:
snowfree_points <- read.csv('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/k_vhf_snowfree.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snowfree_points)
# convert from lat/long to UTM to match rasters
snowfree_points <- st_as_sf(snowfree_points, coords = c("location.long", "location.lat"), crs = 4326)
snowfree_points <- st_transform(snowfree_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snow_points)
plot(st_geometry(snowfree_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

###Calculate peak of quadratic relationship for elevation

summary(sf_glm_8d_us) #get unscaled coefficients for kdem & kdem2
sf_dem_peak = -((8.012*10^-3) / (2 * (-1.697*10^-6)))
sf_dem_peak # peak occurs at 2360.636 metres 
min(sf_train$Band_1)
max(sf_train$Band_1)
summary(sf_train$Band_1) #max value is 1858 but the model is laying optimal elevation is 2360???
summary(sf_test$Band_1) #max value in test dataset is 1951
#Plot Relationship
# Coefficients from the model
beta_0 <- -10.07          # Intercept
beta_1 <- 8.012e-03        # Coefficient for Band_1
beta_2 <- -1.697e-06      # Coefficient for I(Band_1^2)

# Create a sequence of kdem values around the observed range
kdem_values <- seq(500, 4000, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * kdem_values + beta_2 * kdem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(kdem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Late Winter")
abline(v = 2360.64, col = "red", lty = 2)  # Add vertical line at the peak
