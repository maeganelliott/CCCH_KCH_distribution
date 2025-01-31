#script to create seasonal RSF models for clear creek gps data
## EARLY WINTER
#DeMars et al. (2020) used the glmmTMB function from the glmmTMB package
#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r

install.packages("easystats")
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
library(car)
library(AICcmodavg)
library(ggplot2)

set.seed(123)

#Import used and available points with covariates previously extracted to them
ew_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_ew_data_core.csv')

## EARLY WINTER

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(ew_data$landcover2015_250m)
ew_data$landcover2015_250m <- factor(ew_data$landcover2015_250m)

# View the current levels
levels(ew_data$landcover2015_250m)

# Relabel the levels by first creating a vector of the current levels
landcover_levels <- levels(ew_data$landcover2015_250m)
# Modify specific levels by direct replacement
landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10
# Assign the modified levels back to the factor
levels(ew_data$landcover2015_250m) <- landcover_levels
# Verify the new levels
levels(ew_data$landcover2015_250m)

str(ew_data$ccfiredec)
ew_data$ccfiredec <- factor(ew_data$ccfiredec)
levels(ew_data$ccfiredec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(ew_data), 0.8*nrow(ew_data))
ew_train = ew_data[n,]
ew_test = ew_data[-n,] 
#############################################################################
#Univariate environmental glms - unscaled versions - only necessary for covariates where effect of different scales needs to be checked.
#Unique for each season.
#Landcover
ew_lc30_1 <- glm(presence~cclandcover2015, family=binomial, data=ew_train)
ew_lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=ew_train)
ew_lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=ew_train)

summary(ew_lc30_1)
plot(ew_lc30_1)
summary(ew_lc100_1)
plot(ew_lc100_1)
summary(ew_lc250_1)
plot(ew_lc250_1) #LC 250 is the best based on AIC model selection. Low, negative z scores.
#AIC - Land cover
ew_lc_glms <- list(ew_lc30_1, ew_lc100_1, ew_lc250_1)
ew_lc_glms_names <- c('ew_lc30_1', 'ew_lc100_1', 'ew_lc250_1')
aictab(cand.set = ew_lc_glms, modnames = ew_lc_glms_names)

#Lichen
ew_lichbi30 <- glm(presence~cclichen_bi, family=binomial, data=ew_train)
ew_lichenbi100 <- glm(presence~cclichen_bi100, family=binomial, data=ew_train)
ew_lichenbi250 <- glm(presence~cclichen_bi250, family=binomial, data=ew_train)
summary(ew_lich30_1)
summary(ew_lichbi30)# Lichen binomial 30 m is the best based on AIC model selection. High z scores. 
summary(ew_lich100_1)
summary(ew_lich250_1) 
#AIC - lichen PFT
ew_lich_glms <- list(ew_lichbi30, ew_lichenbi100, ew_lichenbi250)
ew_lich_glms_names <- c('ew_lichbi30', 'ew_lichenbi100', 'ew_lichenbi250')
aictab(cand.set = ew_lich_glms, modnames = ew_lich_glms_names)

#Conifer
ew_conbi30 <- glm(presence~ccconifer_bi, family=binomial, data=ew_train)
ew_conbi100 <- glm(presence~ccconifer_bi100, family=binomial, data=ew_train)
ew_conbi250 <- glm(presence~ccconifer_bi250, family=binomial, data=ew_train)
summary(ew_con30_1)
summary(ew_conbi100)# Conifer bivariate 250 m is the best based on AIC model selection. Low z scores. 
summary(ew_con250_1)
#AIC - conifer PFT
ew_con_glms <- list(ew_conbi30, ew_conbi100, ew_conbi250)
ew_con_glms_names <- c('ew_conbi30', 'ew_conbi100', 'ew_conbi250')
aictab(cand.set = ew_con_glms, modnames = ew_con_glms_names)

#Deciduous shrub
ew_decidbi30 <- glm(presence~ccdecid_bi, family=binomial, data=ew_train)
ew_decidbi100 <- glm(presence~ccdecid_bi100, family=binomial, data=ew_train)
ew_decidbi250 <- glm(presence~ccdecid_bi250, family=binomial, data=ew_train)
summary(ew_decidbi30)
summary(ew_decid100_1) #decid bi 30 m is best based on AIC
summary(ew_decid250_1) 
#AIC - deciduous shrub PFT
ew_decid_glms <- list(ew_decidbi30, ew_decidbi100, ew_decidbi250)
ew_decid_glms_names <- c('ew_decidbi30', 'ew_decidbi100', 'ew_decidbi250')
aictab(cand.set = ew_decid_glms, modnames = ew_decid_glms_names)

#graminoid
ew_grambi30 <- glm(presence~ccgram_bi, family=binomial, data=ew_train)
ew_grambi100 <- glm(presence~ccgram_bi100, family=binomial, data=ew_train)
ew_grambi250 <- glm(presence~ccgram_bi250, family=binomial, data=ew_train)
summary(ew_grambi100) #graminoid bivariate 250 m is the best based on AIC model selection. 
ew_gram_glms <- list(ew_grambi30, ew_grambi100, ew_grambi250)
ew_gram_glm_names <- c('ew_grambi30', 'ew_grambi100', 'ew_grambi250')
aictab(cand.set = ew_gram_glms, modnames = ew_gram_glm_names)


#Include: elev (q), slope (q), waterd, landcover250, lichenbi30, conbi250, grambi250, fire

##############################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#https://github.com/rygill/caribou_movement_ecology/tree/main

#Early Winter:
#Full covariate set below: Only include first scale of each repeat variable. 
# Add the transformed terms to the dataset to check collinearity
ew_data$ccslope2 <- ew_data$ccslope^2 #quadratic transformation
ew_data$ccdem2 <- ew_data$ccdem^2

#if the correlation code doesn't work, use this to check for missing columns in dataset.
#colnames(ew_data)
#missing_columns <- setdiff(c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_250m", "cclichen_bi", "ccconifer_bi250", "ccgram_bi250", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_100", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000"), colnames(ew_data))

# Print missing columns
print(missing_columns)

any(c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_250m", "cclichen_bi", "ccconifer_bi250", "ccdecid_bi", "ccgram_bi100", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_100", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000", "road2018_2000") %in% colnames(ew_data))

ew_dat.cor = ew_data[,c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_250m", "cclichen_bi", "ccconifer_bi250", "ccgram_bi250", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")]
names(ew_dat.cor) = c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_250m", "cclichen_bi", "ccconifer_bi250", "ccgram_bi250", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")
cor_matrix <- cor(ew_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(ew_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/cc_distribution/output/stats/ew_correlation_matrix_core.csv", row.names = TRUE)
#slope and roughness highly correlated (0.9)
#Fire, fire decade, and fire sum all highly correlated (0.9)
#distance to disturbance, mining and roads all highly correlated (0.9)
#Buffered roads correlated with buffered disturbance (0.8)
#line density correlated with disturbance density (0.8)
#Conifer percent cover and conifer binomial moderately correlated as are decid percent cover and decid binomial (0.6) but below cutoff. Will still keep the veg bivariate/percent cover variables separate.


####################################################################
#Some experimentation
#core data: elev (q), slope (q), waterd, landcover250, lichenbi30, conbi250, grambi250, fire
#full model
ew_test1 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
ew_test2 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + landcover2015_250m + scale(ccfire2018), family=binomial,data=ew_train)
summary(ew_test2)
#core data:
#In full model only non-vegetated class is not significant and fire has negative non-significant effect
#both grasslands and graminoid have positive effects
#elevation and slope have significant quadratic effects
#Elevation and deciduous forest (negative) are primary drivers of model
#coefficients make sense

#New environmental base model selection:
#include: elev (q), slope (q), waterd, landcover250, lichenbi30, conbi250, grambi250, fire
#exclude waterd - not biologically relevant, layer too coarse.
#excluding deciduous shrub - not biologically relevant for winter
#excluding graminoid - had a small negative coefficient, weak relationship, not contributing very much from a biological perspective.
#exclude fire variables

#all variables
ew_glm_1e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
#all variables minus fire
ew_glm_2e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250), family=binomial,data=ew_train) 
#all variables minus graminoid
ew_glm_3e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
#all variables minus conifer
ew_glm_4e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
#all variables minus lichen
ew_glm_5e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
#all variables minus landcover
ew_glm_6e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train)  
#all variables minus slope
ew_glm_7e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
#all variables minus elevation
ew_glm_8e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train) 
#all variables minus conifer and lichen
#ew_glm_7e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018), family=binomial,data=ew_train)   
#Null
ew_glm_9e <- glm(presence~1, family=binomial,data=ew_train)

ew_e_glms <- list(ew_glm_1e, ew_glm_2e, ew_glm_3e, ew_glm_4e, ew_glm_5e, ew_glm_6e, ew_glm_7e, ew_glm_8e, ew_glm_9e)
ew_e_glms_names <- c('ew_glm_1e', 'ew_glm_2e', 'ew_glm_3e', 'ew_glm_4e', 'ew_glm_5e', 'ew_glm_6e', 'ew_glm_7e', 'ew_glm_8e', 'ew_glm_9e')
aictab(cand.set = ew_e_glms, modnames = ew_e_glms_names)

summary(ew_glm_2e)
#fire is the weakest variable
#full model is second best with delta AIC of 2.01. 
#graminoid is second weakest variable
#elevation and landcover are most important for model fit
#proceed with glm 1e (full mod), including wildfire because of biological relevance in winter season. Very small negative effect.

##################################################
## Model Calibration - environmental base model
#plot predicted vs observed probabilities of occurrence:
#Predict probabilities on the testing data
predicted_probs_test <- predict(ew_glm_2e, newdata = ew_test, type = "response")

#Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = ew_test$presence,       # Binary response variable in the test data
  predicted = predicted_probs_test   # Predicted probabilities from the model
)

#Plot observed vs predicted probabilities
ggplot(results_test, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.5) +    # Scatter plot of predicted vs observed
  geom_smooth(method = "loess", color = "blue") +  # Smoothed trend line
  labs(x = "Predicted Probability", y = "Observed Presence (Binary)", 
       title = "Predicted vs. Observed Probabilities (Testing Data)") +
  theme_minimal()

#Calibration plot
# Step 1: Predict probabilities using the testing data
predicted_probs_test <- predict(ew_glm_2e, newdata = ew_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = ew_test$presence,       # Binary response variable
  predicted = predicted_probs_test   # Predicted probabilities
)

# Step 3: Bin the predicted probabilities into deciles (10 equal-sized bins)
library(dplyr)
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

#Good calibration, the model preforms well until higher predicted probabilities where deviation increases slightly
#Slight over and under prediction across model

#Evaluation with AUC / ROC
#unscaled model:
ew_glm_2e_us <- glm(presence~ccdem + I(ccdem^2) + ccslope + I(ccslope^2) + landcover2015_250m + cclichen_bi + ccconifer_bi250 + ccgram_bi250, family=binomial,data=ew_train) 
#evaluate
ew_glm_2e_eval <- pa_evaluate(predict(ew_glm_2e_us, ew_test[ew_test$presence==1, ]), predict(ew_glm_2e_us, ew_test[ew_test$presence==0, ]))
print(ew_glm_2e_eval)
plot(ew_glm_2e_eval, "ROC")
#AUC=0.765 #Good AUC score, moderate performance
###################
#Model performance with easystats
r2(ew_glm_2e_us) #Tjur's R2 = 0.147

#prediction map:
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m_core.tif') %>%
  resample(ccmask, method="near")
plot(rast250m$landcover2015_250m)

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_250m <- rast250m$landcover2015_250m |> #resample landcover 250m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(18, 16)
cclichen_bi <- envrast$cclichen_bi
ccconifer_bi250 <- rast250m$ccconifer_bi250 |>
  resample(ccmask, method='near')
ccgram_bi250 <- rast250m$ccgram_bi250 |>
  resample(ccmask, method='near')
#ccfire2018 <- firerast$ccfire2018

rasters <- c(ccdem, ccslope, landcover2015_250m, cclichen_bi, ccconifer_bi250, ccgram_bi250) #combine rasters.
st_crs(rasters)
#run predict function on test model
p1 <- predict(rasters, ew_glm_2e_us, type="response")
plot(p1)

#Add gps points to map
#filtered used points w/o covariates extracted to them:
ew_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/GPS_collar_data/cc_gps_2018_earlywinter_core.csv')
st_crs(ew_points)
# convert from lat/long to UTM to match rasters
ew_points <- st_as_sf(ew_points, coords = c("location_long", "location_lat"), crs = 4326)
ew_points <- st_transform(ew_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(ew_points)
plot(st_geometry(ew_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

######################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation
#didn't include larger buffers because of switch in selection based on distributions (dist 2000 & 3000, mine 1000 & 2000 & 3000, road 1000, distden 2000 & 3000, lineden 2000 & 3000)
#disturbance 30 m buffer
ew_glm_1d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(dist2018_30), family=binomial,data=ew_train) 
#disturbance 250 m buffer
ew_glm_2d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(dist2018_250), family=binomial,data=ew_train) 
#disturbance 500 m buffer
ew_glm_3d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(dist2018_500), family=binomial,data=ew_train) 
#disturbance 1000 m buffer 
ew_glm_4d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(dist2018_1000), family=binomial,data=ew_train) 
#dist 1000 m buff * elevation
ew_glm_4di <- glm(presence~scale(ccdem) * scale(dist2018_1000) + scale(I(ccdem^2)) * scale(dist2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250), family=binomial,data=ew_train) 
#mining 30 m buffer
ew_glm_5d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(mine2018_30), family=binomial,data=ew_train) 
#mining 250 m buffer
ew_glm_6d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(mine2018_250), family=binomial,data=ew_train) 
#mining 500 m buffer
ew_glm_7d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(mine2018_500), family=binomial,data=ew_train) 
#road 30 m buffer
ew_glm_8d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(road2018_30), family=binomial,data=ew_train) 
#road 250 m buffer
ew_glm_9d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(road2018_250), family=binomial,data=ew_train) 
#road 500 m buffer
ew_glm_10d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(road2018_500), family=binomial,data=ew_train) 
#road 1000 m buffer
#ew_glm_11d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(road2018_1000), family=binomial,data=ew_train) 
#disturbance density 250 m
ew_glm_12d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(distden2018_250), family=binomial,data=ew_train) 
#disturbance density 500 m
ew_glm_13d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(distden2018_500), family=binomial,data=ew_train) 
#disturbance density 1000 m
ew_glm_14d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(distden2018_1000), family=binomial,data=ew_train) 
#disturbance density 2000 m
#ew_glm_15d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(distden2018_2000), family=binomial,data=ew_train) 
#disturbance density 3000 m
#ew_glm_16d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(distden2018_3000), family=binomial,data=ew_train) 
#line density 250 m
ew_glm_17d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(lineden2018_250), family=binomial,data=ew_train) 
#line density 500 m
ew_glm_18d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(lineden2018_500), family=binomial,data=ew_train) 
#line density 1000 m
ew_glm_19d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(lineden2018_1000), family=binomial,data=ew_train) 
#line density 2000 
#ew_glm_20d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(lineden2018_2000), family=binomial,data=ew_train) 
#line density 2000 * elevation
#ew_glm_20di <- glm(presence~scale(ccdem) * scale(lineden2018_2000) + scale(I(ccdem^2)) * scale(lineden2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250), family=binomial,data=ew_train) 
#line density 3000 m
#ew_glm_21d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(lineden2018_3000), family=binomial,data=ew_train) 
#line density 3000 * elevation
#ew_glm_21di <- glm(presence~scale(ccdem) * scale(lineden2018_3000) + scale(I(ccdem^2)) * scale(lineden2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250), family=binomial,data=ew_train) 
#null
ew_glm_22d <- glm(presence~1,family=binomial,data=ew_train) 
#glm 1e
ew_glm_2e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250), family=binomial, data=ew_train)
#lineden 3000 + road 250 m buff
#ew_glm_21d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250) + scale(ccfire2018) + scale(lineden2018_3000) + scale(road2018_250), family=binomial,data=ew_train) 



ew_d_glms <- list(ew_glm_1d, ew_glm_2d, ew_glm_3d, ew_glm_4d, ew_glm_4di, ew_glm_5d, ew_glm_6d, ew_glm_8d, ew_glm_9d, ew_glm_10d, ew_glm_12d, ew_glm_14d, ew_glm_17d, ew_glm_18d, ew_glm_19d, ew_glm_22d, ew_glm_2e)
ew_d_glms_names <- c('ew_glm_1d', 'ew_glm_2d', 'ew_glm_3d', 'ew_glm_4d', 'ew_glm_4di', 'ew_glm_5d', 'ew_glm_6d', 'ew_glm_8d', 'ew_glm_9d', 'ew_glm_10d', 'ew_glm_12d', 'ew_glm_14d', 'ew_glm_17d', 'ew_glm_18d', 'ew_glm_19d', 'ew_glm_22d', 'ew_glm_2e')
aictab(cand.set = ew_d_glms, modnames = ew_d_glms_names)

#prior to removal due to distribution of used/available points, top variable is line density 3000 m (positive effect), followed by disturbance density 3000 m, then both line and disturbance density 2000 m radii
#dist 1000 m buffer is top buffered variable
#road 500 m buffer is best buffered road variable, but not one of the best models.
#all top disturbance variables are correlated
#The smallest disturbance, road and mining buffers (30 m and 250 m) performed worse than the base environmental model. All have positive but non-significant effects.
#proceed with glm 21d

summary(ew_glm_4di) 
summary(ew_glm_19d)


# Model Calibration
## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(ew_glm_4di, newdata = ew_test, type = "response")

# Step 3: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = ew_test$presence,       # Binary response variable in the test data
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
# Bin the predicted probabilities using deciles.
results_calibration <- results_test %>%
  mutate(predicted_bin = ntile(predicted, 10)) %>%  # Divide into 10 bins
  group_by(predicted_bin) %>%
  summarize(
    mean_predicted = mean(predicted),     # Mean predicted probability in each bin
    observed_rate = mean(observed)        # Proportion of 1s in each bin (observed event rate)
  )

# Create the calibration plot
ggplot(results_calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Perfect calibration line
  labs(x = "Mean Predicted Probability", y = "Observed Event Rate",
       title = "Calibration Plot (Testing Data)") +
  theme_minimal()
#good calibration, the model performs well. Over and under prediction has slightly decreased from base environmental model.

#Test the model / model evaluation
#top model with unscaled predictors
ew_glm_4di_us <- glm(presence~ccdem * dist2018_1000 + I(ccdem^2) * dist2018_1000 + ccslope + I(ccslope^2) + landcover2015_250m + cclichen_bi + ccconifer_bi250 + ccgram_bi250, family=binomial,data=ew_train) 

#evaluate
ew_glm4di_eval <- pa_evaluate(predict(ew_glm_4di_us, ew_test[ew_test$presence==1, ]), predict(ew_glm_4di_us, ew_test[ew_test$presence==0, ]))
print(ew_glm4di_eval)
plot(ew_glm4di_eval, "ROC")
#AUC=0.782 Acceptable/fair discrimination.
###################################################
#Model performance with easystats
r2(ew_glm_4di_us) #Tjur's R2 = 0.166

windows()
check_model(ew_glm_4di_us)
#effect size plot with see pkg from easystats
plot(effectsize(ew_glm_4di_us)) +
  ggplot2::labs(title = "Clear Creek - Early Winter 2018") +
  ggplot2::scale_y_discrete(labels = c(
    "ccdem" = "Elevation",
    "dist2018_1000" = "Disturbance 1000 m",
    "I(ccdem^2)" = "Elevation^2",
    "ccslope" = "Slope",
    "I(ccslope^2)" = "Slope^2",
    "landcover2015_100m5" = "Deciduous/mixed 100 m",
    "landcover2015_100m8" = "Shrublands 100 m",
    "landcover2015_100m10" = "Grasslands 100 m",
    "landcover2015_100m14" = "Wetlands 100 m",
    "landcover2015_100m16" = "Non-vegetated 100 m",
    "cclichen_bi" = "Lichen 30 m",
    "ccconifer_bi250" = "Conifer 250 m",
    "ccgram_bi250" = "Graminoid 250 m",
    "ccdem:dist2018_1000" = "Elevation:Disturbance 1000 m",
    "dist2018_1000:I(ccdem^2)" = "Elevation^2:Disturbance 1000 m"
  )) 

plot(parameters(ew_glm_4di_us)) +
  ggplot2::labs(title = "Early Winter")

#Save models
saveRDS(ew_glm_4di, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_ew_scaled.rds")
saveRDS(ew_glm_4di_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_ew_unscaled.rds")

#Try out predict function 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
#firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccfirerast2018_1.tif')
distbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistbuffrast2018_core.tif')
#distdenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistdenrast2018_1.tif')
#linedenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/cclinedenrast2018_core.tif')
#minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccminebuffrast2018_1.tif')
#roadbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccroadbuffrast2018_1.tif')
#rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m2.tif')%>%
 #resample(ccmask, method="near")
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m_core.tif') %>%
  resample(ccmask, method="near")

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_250m <- rast250m$landcover2015_250m |> #resample landcover 250m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(18, 16)
cclichen_bi <- envrast$cclichen_bi
ccconifer_bi250 <- rast250m$ccconifer_bi250 |>
  resample(ccmask, method='near')
ccgram_bi250 <- rast250m$ccgram_bi250 |>
  resample(ccmask, method='near')
#lineden2018_3000 <- linedenrast$lineden2018_3000
dist2018_1000 <- distbuffrast$dist2018_1000

rasters <- c(ccdem, ccslope, landcover2015_250m, cclichen_bi, ccconifer_bi250, ccgram_bi250, dist2018_1000) #combine rasters. Make sure they are separated by commas.
#run predict function on disturbance model
p1 <- predict(rasters, ew_glm_4di_us, type="response")
plot(p1) #plots!

windows()

#Add gps points to map
#filtered used points w/o covariates extracted to them:
ew_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/GPS_collar_data/cc_gps_2018_earlywinter_core.csv')
st_crs(ew_points)
# convert from lat/long to UTM to match rasters
ew_points <- st_as_sf(ew_points, coords = c("location_long", "location_lat"), crs = 4326)
ew_points <- st_transform(ew_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(ew_points)
plot(st_geometry(ew_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

writeRaster(p1,'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_ew_mod4di_prediction.tif' , overwrite=TRUE)

#Add human disturbances to map
# Read, rasterize, and clip buffered clear creek disturbances 2018.
#30 m buffer
dist2018_30_shp <- st_read('G:/Shared drives/Clear Creek and Klaza Data/Surface Disturbance/Buffered disturbance/CC_dist_2018_30m.shp', quiet=TRUE)
st_crs(dist2018_30_shp)
dist2018_30_utm <- st_as_sf(dist2018_30_shp, coords = c("location_long", "location_lat"), crs = 4326)
dist2018_30_utm <- st_transform(dist2018_30_utm, crs = st_crs(rasters))

plot(st_geometry(dist2018_30_utm), add = TRUE, col = "red", pch = 19, cex = 0.6) 

writeRaster(p1,'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_ew_mod4di_predict_pts_dist.tif' , overwrite=TRUE)

#####################################################################################

####Calculate peak of quadratic relationship for elevation
windows()
summary(ew_glm_4di_us) #get unscaled coefficients for elevation
#dem = 0.02972, dem2 = -1.239*10^-5
ew_dem_peak = -(0.02972 / (2 * (-1.239*10^-5)))
ew_dem_peak # peak occurs at 1199.354 metres 

#Plot Relationship
# Coefficients from the model
beta_0 <- -21.23          # Intercept
beta_1 <- 0.02972         # Coefficient for dem
beta_2 <- -1.239e-05      # Coefficient for I(dem^2)

# Create a sequence of dem values around the observed range
dem_values <- seq(500, 2000, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * dem_values + beta_2 * dem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(dem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Early Winter")
abline(v = 1199.354, col = "red", lty = 2)  # Add vertical line at the peak

#########################################################
#Model Figures
#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r
install.packages("sjPlot")
library(sjPlot)
install.packages("gridExtra")
library(gridExtra)
#Model ew_glm_4di
ew_glm_4di <- glm(presence~scale(ccdem) * scale(dist2018_1000) + scale(I(ccdem^2)) * scale(dist2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_250m + scale(cclichen_bi) + scale(ccconifer_bi250) + scale(ccgram_bi250), family=binomial,data=ew_train) 

sjPlot::plot_model(ew_glm_4di,
           type = "pred",
           show.data = T,
           terms = c("ccdem", "ccdem^2"))


ggplot(data = ew_train + 
  geom_boxplot(aes(x = ccdem, y = presence)))

#Generate a figure of the modelling results
blue1 = "#4774A0"
blue2 = "#7BA2C9"
blueA = "#BCD5EE"
red = "#FF0000"

ggplot() +
  geom_point(data = dat.hr, aes(x = herd, y = HR, color = period)) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("During Lockdown", "1 Year Prior", "2 Years Prior", "After"),
                      name = "period")


#get predicted values for panel e to order the axis:
per.table = get_model_data(sel_mod, type = c('pred'))[[1]]
per.table$period = as.factor(ifelse(per.table$x == 1, '2021',
                                    ifelse(per.table$x == 2, '2020',
                                           ifelse(per.table$x == 3, '2019','2022'))))
per.table = per.table[,c(8,2,4,5)]

#Generate plots
c <- 
  plot_model(ew_glm_4di,
             type = "pred",
             line.size = 0.2,
             terms = c("Elevation")) +
  geom_point(data = dat.hr, aes(x = ccdem, y = presence, col = period), alpha = 0.6, pch = 16, size = 0.5) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("2021", "2020", "2019", "2022"),
                      name = "period") + 
  ggtitle("c)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  #scale_x_continuous(limits = c(0,1), expand = c(0,0.02)) + 
  #scale_y_continuous(limits = c(0,2100), expand = c(0,0.1)) + 
  ylab(expression(paste("Predicted home range area (",km^2,")"))) +
  xlab("Standardised slope") +
  coord_cartesian(ylim =  c(0,2100))