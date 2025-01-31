#VHF GLM model selection
#Oct 23, 2024
#Reviewed Nov 21, 2024

library(easystats)
library(see)
library(effectsize)
library(parameters)
library(DHARMa)
library(patchwork)
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

##SNOW SEASON

#Import used and available points with covariates previously extracted to them
snow_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_snow_data.csv')
#check data:
str(snow_data) #numeric and integers, will work with scale function
snow_data$Band_1 <- as.numeric(snow_data$Band_1) #set elevation as numeric
colSums(is.na(snow_data)) #check for NAs - 1 point with NAs to remove. 
snow_data <- na.omit(snow_data)


#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(snow_data), 0.8*nrow(snow_data))
snow_train = snow_data[n,]
snow_test = snow_data[-n,] 

colSums(is.na(snow_train))
colSums(is.na(snow_test))

#######################################################################
#Run univariate models to determine best scale for each PFT veg raster
#Lichen
lichbi30 <- glm(presence~lichen_bi00, family=binomial, data=snow_train)
lichenbi100 <- glm(presence~lichen_bi00_100, family=binomial, data=snow_train)
lichenbi250 <- glm(presence~lichen_bi00_250, family=binomial, data=snow_train)
summary(lichenbi250) # Lichen bivariate 30 m is the best based on AIC model selection. Medium, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi00, family=binomial, data=snow_train)
conbi100 <- glm(presence~conifer_bi00_100, family=binomial, data=snow_train)
conbi250 <- glm(presence~conifer_bi00_250, family=binomial, data=snow_train)
summary(conbi100) # Conifer bivariate 100 m is the best based on AIC model selection but all scales very similar. Positive effect.  
#AIC 
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
#decidbi30 <- glm(presence~decid_bi00, family=binomial, data=snow_train)
#decidbi100 <- glm(presence~decid_bi00_100, family=binomial, data=snow_train)
#decidbi250 <- glm(presence~decid_bi00_250, family=binomial, data=snow_train)
#summary(decidbi100) # decid bivariate 100 m is the best based on AIC model selection but all scales very similar. Positive effect.  
#AIC 
#decid_glms <- list(decidbi30, decidbi100, decidbi250)
#decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
#aictab(cand.set = decid_glms, modnames = decid_glms_names)

#Graminoid
grambi30 <- glm(presence~gram_bi00, family=binomial, data=snow_train)
grambi100 <- glm(presence~gram_bi00_100, family=binomial, data=snow_train)
grambi250 <- glm(presence~gram_bi00_250, family=binomial, data=snow_train)
summary(grambi100) # Graminoid bivariate 100 m is the best based on AIC model selection but all scales very similar. Negative effect.  
#AIC 
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glms_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glms_names)

#Include: dem, slope, lichenbi30, coniferbi100, grambi100, fire
#exclude deciduous shrub for winter?

#Check Correlation of Snow season variables
#add quadratic terms
snow_data$ccslope2 <- snow_data$ccslope^2 #quadratic transformation
snow_data$ccdem2 <- snow_data$Band_1^2 #remember elevation is called Band_1

#missing_cols <- setdiff(c("Band_1", "ccdem2", "ccslope", "ccslope2", "lichen_bi00", "lichen_bi00_100", "lichen_bi00_250", "conifer_bi00", "conifer_bi00_100", "conifer_bi00_250", "decid_bi00_100", "gram_bi00_100", "fire2000", "dist2000_30", "dist2000_500", "dist2000_1000", "dist2000_2000", "dist2000_3000", "lineden2000_250", "lineden2000_500", "lineden2000_1000", "lineden2000_2000", "lineden2000_3000"), colnames(snow_data))
#missing_cols # check for miss-matched names if correlation code doesn't work.


snow_dat.cor = snow_data[,c("Band_1", "ccdem2", "ccslope", "ccslope2", "lichen_bi00", "lichen_bi00_100", "lichen_bi00_250", "conifer_bi00", "conifer_bi00_100", "conifer_bi00_250", "decid_bi00_100", "gram_bi00_100", "fire2000", "dist2000_30", "dist2000_500", "dist2000_1000", "dist2000_2000", "dist2000_3000", "lineden2000_250", "lineden2000_500", "lineden2000_1000", "lineden2000_2000", "lineden2000_3000")]
names(snow_dat.cor) = c("Band_1", "ccdem2", "ccslope", "ccslope2", "lichen_bi00", "lichen_bi00_100", "lichen_bi00_250", "conifer_bi00", "conifer_bi00_100", "conifer_bi00_250", "decid_bi00_100", "gram_bi00_100", "fire2000", "dist2000_30", "dist2000_500", "dist2000_1000", "dist2000_2000", "dist2000_3000", "lineden2000_250", "lineden2000_500", "lineden2000_1000", "lineden2000_2000", "lineden2000_3000")
cor_matrix <- cor(snow_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(snow_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/cc_distribution/output/stats/vhf_snow_correlation_matrix2.csv", row.names = TRUE)
#some correlation between disturbance variables to watch for. 
#conifer and elevation negatively correlated

#######################################################################

##Snow Season GLMs

#Some experimentation. Note elevation and conifer are negatively correlated. 
#variables: dem, slope, lichenbi250, coniferbi250, grambi2100, fire
#test with elevation
# all variables - dem, slope, lichenbi250 grambi100, fire (except for conifer)
s_test1 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi00_250) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#remove gram
s_test2 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi00_250) + scale(fire2000), family=binomial,data=snow_train)
#remove gram and fire
s_test3 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi_00250), family=binomial,data=snow_train)
#elevation, lichen, fire
s_test4 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(lichen_bi00_250) + scale(fire2000), family=binomial,data=snow_train)
#all variables with linear slope term
s_test5 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(lichen_bi00_250) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#fire only
s_test6 <- glm(presence~scale(fire2000), family=binomial, data=snow_train)
#fire and slope
s_test7 <- glm(presence~scale(ccslope) + scale(fire2000), family=binomial, data=snow_train)
#slope only
s_test8 <- glm(presence~scale(ccslope), family=binomial,data=snow_train)
#test with conifer
#all variables (minus elevation)
s_test9 <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#lichen, conifer and fire (no slope)
s_test10 <- glm(presence~scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#lichen and conifer
s_test11 <- glm(presence~scale(lichen_bi00_250) + scale(conifer_bi00_100), family=binomial,data=snow_train)
#lichen, conifer, slope
s_test12 <- glm(presence~scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(ccslope), family=binomial,data=snow_train)

summary(s_test9)

#elevation has the largest effect size, significant quadratic relationship
#slope has non-significant quadratic relationship. When linear term modelled on its own it has a significant negative effect. Proceed with linear slope term.
#lichen has positive non-significant effect in full mod (with elevation)
#graminoid has negative non-significant effect in full mod (with elevation). Removing graminoid brings lichen closer to significance.
#fire has smallest effect size, positive and non-significant when modelled with elevation. Negative and non-significant when modelled with conifer. 
#fire effect direction is negative when modelled with lichen 30 m and positive when modelled with lichen 250 m. 
#lichen 250 m is significant and lichen 30 m is not significant (positive effects). 
#conifer has a positive significant effect
#need to use AIC to choose between elevation and conifer

#GLM - based on covariate exploration / distributions
#Include: dem, slope (linear), lichenbi30, coniferbi250, grambi250, fire

#conifer mods:
#global with conifer
snow_glm_1e <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#all variables minus fire
snow_glm_2e  <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100), family=binomial,data=snow_train)
#all variables minus graminoid
snow_glm_3e <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#all variables minus slope
snow_glm_4e <- glm(presence~scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#all variables minus slope and graminoid
#snow_glm_5e <- glm(presence~scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#elevation mods:
#global with elevation
snow_glm_5e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(lichen_bi00_250) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#primary model drivers (lichen and elevation)
#snow_glm_8e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(lichen_bi00_250), family=binomial,data=snow_train)
#all variables minus fire
snow_glm_6e  <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(lichen_bi00_250) + scale(gram_bi00_100), family=binomial,data=snow_train)
#all variables minus graminoid
snow_glm_7e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(lichen_bi00_250) + scale(fire2000), family=binomial,data=snow_train)
#all variables minus slope
snow_glm_8e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(lichen_bi00_250) + scale(gram_bi00_100) + scale(fire2000), family=binomial,data=snow_train)
#all variables minus slope and graminoid
#snow_glm_10e <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(lichen_bi00_250) + scale(fire2000), family=binomial,data=snow_train)
#null
snow_glm_9e <- glm(presence~1, family=binomial,data=snow_train)

#AIC Model Selection
snow_e_glms <- list(snow_glm_1e, snow_glm_2e, snow_glm_3e, snow_glm_4e, snow_glm_5e, snow_glm_6e, snow_glm_7e, snow_glm_8e, snow_glm_9e)
snow_e_glms_names <- c('snow_glm_1e', 'snow_glm_2e', 'snow_glm_3e', 'snow_glm_4e', 'snow_glm_5e', 'snow_glm_6e', 'snow_glm_7e', 'snow_glm_8e', 'snow_glm_9e')
aictab(cand.set = snow_e_glms, modnames = snow_e_glms_names)

summary(snow_glm_2e)


## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(snow_glm_2e, newdata = snow_test, type = "response")

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
predicted_probs_test <- predict(snow_glm_2e, newdata = snow_test, type = "response")

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

#model tends to over and underpredict due to small sample size. Lots of deviation.

#Evaluation with AUC / ROC
#unscaled model:
snow_glm_2e_us <- glm(presence~ccslope + lichen_bi00_250 + conifer_bi00_100 + gram_bi00_100, family=binomial,data=snow_train)
#evaluate
snow_glm_2e_eval <- pa_evaluate(predict(snow_glm_2e_us, snow_test[snow_test$presence==1, ]), predict(snow_glm_2e_us, snow_test[snow_test$presence==0, ]))
print(snow_glm_2e_eval)
plot(snow_glm_2e_eval, "ROC")
#AUC=0.58 #poor discrimination 
#######################################
#Model performance with easystats
r2(snow_glm_2e_us) #Tjur's R2 = 0.055


#Caribou occurrence prediction map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/cc_vhf_mask_30m_3.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_30m_2.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_100m_2.tif') %>%
  resample(mask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_250m_2.tif') %>%
  resample(mask, method="near")

#prepare rasters
ccslope <- rast30m$ccslope
lichen_bi00_250 <- rast250m$lichen_bi00_250 |>
  resample(mask, method='near')
conifer_bi00_100 <- rast100m$conifer_bi00_100 |>
  resample(mask, method='near')
gram_bi00_100 <- rast100m$gram_bi00_100 |>
  resample(mask, method='near')

rasters <- c(ccslope, lichen_bi00_250, conifer_bi00_100, gram_bi00_100) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, snow_glm_2e_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_snow_mod2e_prediction.tif')

#Add vhf points to map
#filtered used points w/o covariates extracted to them:
snow_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/VHF_collar_data/cc_vhf_snow_cleaned2.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snow_points)
# convert from lat/long to UTM to match rasters
snow_points <- st_as_sf(snow_points, coords = c("location_long", "location_lat"), crs = 4326)
snow_points <- st_transform(snow_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snow_points)
plot(st_geometry(snow_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 


## Phase 2 - Disturbance Model Selection
#double check collinearity before running combined models
#Buffered disturbance 30 m
snow_glm_1d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(dist2000_30), family=binomial,data=snow_train)
#buffered disturbance 500 m
snow_glm_2d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(dist2000_500), family=binomial,data=snow_train)
#buffered disturbance 1000 m
snow_glm_3d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(dist2000_1000), family=binomial,data=snow_train)
#buffered disturbance 2000 m
snow_glm_4d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(dist2000_2000), family=binomial,data=snow_train)
#buffered disturbance 3000 m
snow_glm_5d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(dist2000_3000), family=binomial,data=snow_train)
#line density 250 
snow_glm_6d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(lineden2000_250), family=binomial,data=snow_train)
#line density 500 
snow_glm_7d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(lineden2000_500), family=binomial,data=snow_train)
#line density 1000
snow_glm_8d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(lineden2000_1000), family=binomial,data=snow_train)
#line density 2000
snow_glm_9d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(lineden2000_2000), family=binomial,data=snow_train)
#line density 3000
snow_glm_10d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(lineden2000_3000), family=binomial,data=snow_train)
#buffered dist 3000 + lineden 1000
snow_glm_11d <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100) + scale(dist2000_3000) + scale(lineden2000_1000), family=binomial,data=snow_train)
#glm 2e
snow_glm_2e <- glm(presence~scale(ccslope) + scale(lichen_bi00_250) + scale(conifer_bi00_100) + scale(gram_bi00_100), family=binomial,data=snow_train)
#null
snow_glm_12d <- glm(presence~1, family=binomial,data=snow_train)

#AIC Model Selection
snow_d_glms <- list(snow_glm_1d, snow_glm_2d, snow_glm_3d, snow_glm_4d, snow_glm_5d, snow_glm_6d, snow_glm_7d, snow_glm_8d, snow_glm_9d, snow_glm_10d, snow_glm_2e, snow_glm_11d, snow_glm_12d)
snow_d_glms_names <- c('snow_glm_1d', 'snow_glm_2d', 'snow_glm_3d', 'snow_glm_4d', 'snow_glm_5d', 'snow_glm_6d', 'snow_glm_7d', 'snow_glm_8d', 'snow_glm_9d', 'snow_glm_10d', 'snow_glm_2e', 'snow_glm_11d', 'snow_glm_12d')
aictab(cand.set = snow_d_glms, modnames = snow_d_glms_names)

summary(snow_glm_5d)

#larger buffers have a significant positive effect and improved model fit the most
#smaller buffers have non-significant positive effects
#linear feature density has non-significant positive effect
#tried combining the best performing linear feature density variable with the best performing buffered disturbance variable, but performed worse than the model with only the buffered disturbance variable
#buffered disturbance (3000 m) performed the best. 
#some of the linear feature density variables performed worse than the base environmental model (250, 500, 2000 m radii)
#proceed with glm 5d, with buffered disturbance (3000 m)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(snow_glm_5d, newdata = snow_test, type = "response")

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
predicted_probs_test <- predict(snow_glm_5d, newdata = snow_test, type = "response")

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

#model tends to over and underpredict due to small sample size. Lots of deviation.

#Evaluation with AUC / ROC
#unscaled model:
snow_glm_5d_us <- glm(presence~ccslope + lichen_bi00_250 + conifer_bi00_100 + gram_bi00_100 + dist2000_3000, family=binomial,data=snow_train)
#evaluate
snow_glm_5d_eval <- pa_evaluate(predict(snow_glm_5d_us, snow_test[snow_test$presence==1, ]), predict(snow_glm_5d_us, snow_test[snow_test$presence==0, ]))
print(snow_glm_5d_eval)
plot(snow_glm_5d_eval, "ROC")
#AUC=0.647 #moderate discrimination 
#######################################
#Model performance with easystats
r2(snow_glm_5d_us) #Tjur's R2 = 0.074

windows()
check_model(snow_glm_5d_us)
#effect size plot with see pkg from easystats
plot(effectsize(snow_glm_5d_us)) +
  ggplot2::labs(title = "Clear Creek - Snow 2000") +
  ggplot2::scale_y_discrete(labels = c(
    "ccslope" = "Slope",
    "lichen_bi00_250" = "Lichen 250 m",
    "conifer_bi00_100" = "Conifer 100 m",
    "gram_bi00_100" = "Graminoid 100 m",
    "dist2000_3000" = "Disturbance 3000 m"
  )) 


plot(parameters(snow_glm_5d_us)) +
  ggplot2::labs(title = "Snow")

#Save models
saveRDS(snow_glm_5d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_snow_scaled.rds")
saveRDS(snow_glm_5d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_snow_unscaled.rds")

#Caribou occurrence prediction map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/cc_vhf_mask_30m_3.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_30m_2.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_100m_2.tif') %>%
  resample(mask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_250m_2.tif') %>%
  resample(mask, method="near")
distrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistbuffrast2000_2.tif')
linerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/cclinedenrast2000_2.tif')

#prepare rasters
ccslope <- rast30m$ccslope
lichen_bi00_250 <- rast250m$lichen_bi00_250 |>
  resample(mask, method='near')
conifer_bi00_100 <- rast100m$conifer_bi00_100 |>
  resample(mask, method='near')
gram_bi00_100 <- rast100m$gram_bi00_100 |>
  resample(mask, method='near')
dist2000_3000 <- distrast$dist2000_3000
#ineden2000_1000 <- linerast$lineden2000_1000

rasters <- c(ccslope, lichen_bi00_250, conifer_bi00_100, gram_bi00_100, dist2000_3000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, snow_glm_5d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_snow_mod5d_prediction.tif')

#Add vhf points to map
#filtered used points w/o covariates extracted to them:
snow_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/VHF_collar_data/cc_vhf_snow_cleaned2.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snow_points)
# convert from lat/long to UTM to match rasters
snow_points <- st_as_sf(snow_points, coords = c("location_long", "location_lat"), crs = 4326)
snow_points <- st_transform(snow_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snow_points)
plot(st_geometry(snow_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

#######################################################################################
##SNOW FREE SEASON

#Import used and available points with covariates previously extracted to them
snowfree_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_snowfree_data_2.csv')
#check data:
str(snowfree_data) #numeric and integers, will work with scale function
snowfree_data$Band_1 <- as.numeric(snowfree_data$Band_1) #set elevation as numeric
colSums(is.na(snowfree_data)) #check for NAs - one
snowfree_data <- na.omit(snowfree_data)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(snowfree_data), 0.8*nrow(snowfree_data))
sf_train = snowfree_data[n,]
sf_test = snowfree_data[-n,] 

########################################################################
#Run univariate models to determine best scale for each PFT veg raster
#Lichen
lichbi30 <- glm(presence~lichen_bi00, family=binomial, data=sf_train)
lichenbi100 <- glm(presence~lichen_bi00_100, family=binomial, data=sf_train)
lichenbi250 <- glm(presence~lichen_bi00_250, family=binomial, data=sf_train)
summary(lichenbi250) # Lichen bivariate 30 m is the best based on AIC model selection. Medium, positive z score. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi00, family=binomial, data=sf_train)
conbi100 <- glm(presence~conifer_bi00_100, family=binomial, data=sf_train)
conbi250 <- glm(presence~conifer_bi00_250, family=binomial, data=sf_train)
summary(conbi100) # Conifer bivariate 100 m is the best based on AIC model selection but all scales very similar. Negative effect.  
#AIC 
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~decid_bi00, family=binomial, data=sf_train)
decidbi100 <- glm(presence~decid_bi00_100, family=binomial, data=sf_train)
decidbi250 <- glm(presence~decid_bi00_250, family=binomial, data=sf_train)
summary(decidbi100) # decid bivariate 30 m is the best based on AIC model selection but all scales very similar. Negative  effect.  
#AIC 
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#Graminoid
grambi30 <- glm(presence~gram_bi00, family=binomial, data=sf_train)
grambi100 <- glm(presence~gram_bi00_100, family=binomial, data=sf_train)
grambi250 <- glm(presence~gram_bi00_250, family=binomial, data=sf_train)
summary(grambi30) # Graminoid bivariate 30 m is the best based on AIC model selection but all scales very similar & not significant. Positive effect.  
#AIC 
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glms_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glms_names)

#Include: dem, slope, lichenbi30, coniferbi100, decidbi30, grambi30
#conifer and elevation are negatively correlated
#Check Correlation of Snow season variables
#add quadratic terms
snowfree_data$ccslope2 <- snowfree_data$ccslope^2 #quadratic transformation
snowfree_data$ccdem2 <- snowfree_data$Band_1^2 #remember elevation is called Band_1

#missing_cols <- setdiff(c("Band_1", "ccdem2", "ccslope", "ccslope2", "lichen_bi00", "conifer_bi00_100", "decid_bi00_100", "gram_bi00", "fire2000", "dist2000_30", "dist2000_500", "dist2000_1000", "dist2000_2000", "dist2000_3000", "lineden2000_250", "lineden2000_500", "lineden2000_1000", "lineden2000_2000", "lineden2000_3000"), colnames(snowfree_data))
#missing_cols # check for miss-matched names if correlation code doesn't work.


snowfree_dat.cor = snowfree_data[,c("Band_1", "ccdem2", "ccslope", "ccslope2", "lichen_bi00", "conifer_bi00", "conifer_bi00_100", "conifer_bi00_250", "decid_bi00", "gram_bi00", "fire2000", "dist2000_30", "dist2000_500", "dist2000_1000", "dist2000_2000", "dist2000_3000", "lineden2000_250", "lineden2000_500", "lineden2000_1000", "lineden2000_2000", "lineden2000_3000")]
names(snowfree_dat.cor) = c("Band_1", "ccdem2", "ccslope", "ccslope2", "lichen_bi00", "conifer_bi00", "conifer_bi00_100", "conifer_bi00_250", "decid_bi00", "gram_bi00", "fire2000", "dist2000_30", "dist2000_500", "dist2000_1000", "dist2000_2000", "dist2000_3000", "lineden2000_250", "lineden2000_500", "lineden2000_1000", "lineden2000_2000", "lineden2000_3000")
cor_matrix <- cor(snowfree_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(snow_datfree.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/cc_distribution/output/stats/vhf_snowfree_correlation_matrix.csv", row.names = TRUE)
#some correlation between disturbance variables to watch for. 

#####################################################################
#Snow Free GLMs

#variables: dem, slope, lichenbi30, coniferbi100, decidbi30, grambi30
#elevation and conifer are negatively correlated

#Some experimentation
#all variables with elevation: 
sf_test1 <- glm(presence~scale(Band_1) + scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi00) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#all variables with conifer:
sf_test2 <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#remove decid
sf_test3 <- glm(presence~scale(Band_1) + scale(I(Band_1^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(lichen_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#elevation only with quadratic term
sf_test4 <-glm(presence~scale(Band_1) + scale(I(Band_1^2)), family=binomial, data=sf_train)
#elevation only, linear term.
sf_test5 <-glm(presence~scale(Band_1), family=binomial, data=sf_train)
#slope only
sf_test6 <- glm(presence~scale(ccslope) + scale(I(ccslope^2)), family=binomial, data=sf_train)


summary(sf_test1)

#in full models lichen has a significant positive effect, graminoid has significant positive effect, linear slope term has significant positive effect
#deciduous shrub is weakest predictor
#conifer has negative but non-significant effect in combined model
#elevation has non-significant u-shaped relationship in combined model
#in elevation only model neither elevation term significant, both are positive. 
#when only the linear elevation term is modeled it is very significant and has a positive effect
#In slope only model only the linear term is significant. 
#use AIC model selection to choose between elevation and conifer, and to decide wether to include slope and deciduous shrub in the final model. 

#Environmental model selection
#include: dem, slope, lichenbi30, coniferbi100, decidbi30, grambi30

#Conifer mods
#global 
sf_glm_1e <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#no gram
sf_glm_2e <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(decid_bi00), family=binomial,data=sf_train)
#no decid
sf_glm_3e <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00), family=binomial,data=sf_train)
#no conifer or elevation
sf_glm_4e <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#no lichen
sf_glm_5e <- glm(presence~scale(ccslope) + scale(conifer_bi00_100) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#no slope
sf_glm_6e <- glm(presence~scale(lichen_bi00) + scale(conifer_bi00_100) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#elevation mods
#global
sf_glm_7e <- glm(presence~scale(Band_1) + scale(ccslope) + scale(lichen_bi00) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#no gram
sf_glm_8e <- glm(presence~scale(Band_1) + scale(ccslope) + scale(lichen_bi00) + scale(decid_bi00), family=binomial,data=sf_train)
#no decid
sf_glm_9e <- glm(presence~scale(Band_1) + scale(ccslope) + scale(lichen_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#no lichen
sf_glm_10e <- glm(presence~scale(Band_1) + scale(ccslope) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#no slope
sf_glm_11e <- glm(presence~scale(Band_1) + scale(lichen_bi00) + scale(decid_bi00) + scale(gram_bi00), family=binomial,data=sf_train)
#null
sf_glm_12e <- glm(presence~1, family=binomial, data=sf_train)

#AIC Model Selection
sf_e_glms <- list(sf_glm_1e, sf_glm_2e, sf_glm_3e, sf_glm_4e, sf_glm_5e, sf_glm_6e, sf_glm_7e, sf_glm_8e, sf_glm_9e, sf_glm_10e, sf_glm_11e, sf_glm_12e)
sf_e_glms_names <- c('sf_glm_1e', 'sf_glm_2e', 'sf_glm_3e', 'sf_glm_4e', 'sf_glm_5e', 'sf_glm_6e', 'sf_glm_7e', 'sf_glm_8e', 'sf_glm_9e', 'sf_glm_10e', 'sf_glm_11e', 'sf_glm_12e')
aictab(cand.set = sf_e_glms, modnames = sf_e_glms_names)

summary(sf_glm_7e)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(sf_glm_3e, newdata = sf_test, type = "response")

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
predicted_probs_test <- predict(sf_glm_3e, newdata = sf_test, type = "response")

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

#model tends to over and underpredict due to small sample size. Lots of deviation.

#Evaluation with AUC / ROC
#unscaled model:
sf_glm_3e_us <- glm(presence~ccslope + lichen_bi00 + conifer_bi00_100 + gram_bi00, family=binomial,data=sf_train)
#evaluate
sf_glm_3e_eval <- pa_evaluate(predict(sf_glm_3e_us, sf_test[sf_test$presence==1, ]), predict(sf_glm_3e_us, sf_test[sf_test$presence==0, ]))
print(sf_glm_3e_eval)
plot(sf_glm_3e_eval, "ROC")
#AUC=0.529 #poor discrimination 
###########################################
#Model performance with easystats
r2(sf_glm_3e_us) #Tjur's R2 = 0.054

#Caribou occurrence prediction map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/cc_vhf_mask_30m_3.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_30m_2.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_100m_2.tif') %>%
  resample(mask, method="near")

#prepare rasters
ccslope <- rast30m$ccslope
lichen_bi00 <- rast30m$lichen_bi00
conifer_bi00_100 <- rast100m$conifer_bi00_100 |>
  resample(mask, method='near')
gram_bi00 <- rast30m$gram_bi00

rasters <- c(ccslope, lichen_bi00, conifer_bi00_100, gram_bi00) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, sf_glm_3e_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_snow_mod5d_prediction.tif')

#Add vhf points to map
#filtered used points w/o covariates extracted to them:
snowfree_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/VHF_collar_data/cc_vhf_snowfree_cleaned2.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snowfree_points)
# convert from lat/long to UTM to match rasters
snowfree_points <- st_as_sf(snowfree_points, coords = c("location_long", "location_lat"), crs = 4326)
snowfree_points <- st_transform(snowfree_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snowfree_points)
plot(st_geometry(snowfree_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

##Disturbance model selection
#check correlation before running combined models
#disturbance 30 m buff
sf_glm_1d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(dist2000_30), family=binomial,data=sf_train)
#disturbance 500 m buff
sf_glm_2d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(dist2000_500), family=binomial,data=sf_train)
#disturbance 1000 m buff
sf_glm_3d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(dist2000_1000), family=binomial,data=sf_train)
#disturbance 2000 m buff
sf_glm_4d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(dist2000_2000), family=binomial,data=sf_train)
#disturbance 3000 m buff
sf_glm_5d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(dist2000_3000), family=binomial,data=sf_train)
#lineden 250 m 
sf_glm_6d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(lineden2000_250), family=binomial,data=sf_train)
#lineden 500 m
sf_glm_7d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(lineden2000_500), family=binomial,data=sf_train)
#lineden 1000 m
sf_glm_8d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(lineden2000_1000), family=binomial,data=sf_train)
#lineden 2000 m
sf_glm_9d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(lineden2000_2000), family=binomial,data=sf_train)
#lineden 3000 m
sf_glm_10d <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00) + scale(lineden2000_3000), family=binomial,data=sf_train)
#glm 3e
sf_glm_3e <- glm(presence~scale(ccslope) + scale(lichen_bi00) + scale(conifer_bi00_100) + scale(gram_bi00), family=binomial,data=sf_train)
#null
sf_glm_11d <- glm(presence~1, family=binomial,data=sf_train)

#AIC Model Selection
sf_d_glms <- list(sf_glm_1d, sf_glm_2d, sf_glm_3d, sf_glm_4d, sf_glm_5d, sf_glm_6d, sf_glm_7d, sf_glm_8d, sf_glm_9d, sf_glm_10d, sf_glm_11d, sf_glm_3e)
sf_d_glms_names <- c('sf_glm_1d', 'sf_glm_2d', 'sf_glm_3d', 'sf_glm_4d', 'sf_glm_5d', 'sf_glm_6d', 'sf_glm_7d', 'sf_glm_8d', 'sf_glm_9d', 'sf_glm_10d', 'sf_glm_11d', 'sf_glm_3e')
aictab(cand.set = sf_d_glms, modnames = sf_d_glms_names)

summary(sf_glm_9d)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Predict probabilities on the testing data
predicted_probs_test <- predict(sf_glm_9d, newdata = sf_test, type = "response")

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
predicted_probs_test <- predict(sf_glm_9d, newdata = sf_test, type = "response")

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

#model tends to over and underpredict due to small sample size. Lots of deviation.

#Evaluation with AUC / ROC
#unscaled model:
sf_glm_9d_us <- glm(presence~ccslope + lichen_bi00 + conifer_bi00_100 + gram_bi00 + lineden2000_2000, family=binomial,data=sf_train)
#evaluate
sf_glm_9d_eval <- pa_evaluate(predict(sf_glm_9d_us, sf_test[sf_test$presence==1, ]), predict(sf_glm_9d_us, sf_test[sf_test$presence==0, ]))
print(sf_glm_9d_eval)
plot(sf_glm_9d_eval, "ROC")
#AUC=0.537 #poor discrimination, slight improvement from base model.  
############################################
#Model performance with easystats
r2(sf_glm_9d_us) #Tjur's R2 = 0.069

windows()
check_model(sf_glm_9d_us)
#effect size plot with see pkg from easystats
plot(effectsize(sf_glm_9d_us)) +
  ggplot2::labs(title = "Clear Creek - Snow-free 2000") +
  ggplot2::scale_y_discrete(labels = c(
    "ccslope" = "Slope",
    "lichen_bi00" = "Lichen 30 m",
    "conifer_bi00_100" = "Conifer 100 m",
    "gram_bi00" = "Graminoid 30 m",
    "lineden2000_2000" = "Line density 2000 m"
  )) 

plot(parameters(sf_glm_9d_us)) +
  ggplot2::labs(title = "Snow-free")

#Save models
saveRDS(sf_glm_9d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_snowfree_scaled.rds")
saveRDS(sf_glm_9d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_snowfree_unscaled.rds")

#Caribou occurrence prediction map
#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/vhf_rasters/cc_vhf_mask_30m_3.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
rast30m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_30m_2.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2000_100m_2.tif') %>%
  resample(mask, method="near")
linerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/cclinedenrast2000_2.tif')

#prepare rasters
ccslope <- rast30m$ccslope
lichen_bi00 <- rast30m$lichen_bi00
conifer_bi00_100 <- rast100m$conifer_bi00_100 |>
  resample(mask, method='near')
gram_bi00 <- rast30m$gram_bi00
lineden2000_2000 <- linerast$lineden2000_2000

rasters <- c(ccslope, lichen_bi00, conifer_bi00_100, gram_bi00, lineden2000_2000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, sf_glm_9d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_snowfree_mod9d_prediction.tif')

#Add vhf points to map
#filtered used points w/o covariates extracted to them:
snowfree_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/VHF_collar_data/cc_vhf_snowfree_cleaned2.csv') %>%
  st_as_sf(coords=c('location_long', 'location_lat'), crs=4326) %>%
  st_transform(3578)
st_crs(snowfree_points)
# convert from lat/long to UTM to match rasters
snowfree_points <- st_as_sf(snowfree_points, coords = c("location_long", "location_lat"), crs = 4326)
snowfree_points <- st_transform(snowfree_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(snowfree_points)
plot(st_geometry(snowfree_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 
