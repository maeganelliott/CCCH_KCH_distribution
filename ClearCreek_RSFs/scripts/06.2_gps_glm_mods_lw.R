#script to create seasonal RSF models for clear creek gps data
## LATE WINTER
#DeMars et al. (2020) used the glmmTMB function from the glmmTMB package
#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r

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
library(dplyr)
library(ggplot2)
library(car)
library(AICcmodavg)

set.seed(123)

#Import used and available points with covariates previously extracted to them
lw_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_lw_data_core.csv')

## LATE WINTER

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(lw_data$landcover2015_100m)
lw_data$landcover2015_100m <- factor(lw_data$landcover2015_100m)
levels(lw_data$landcover2015_100m) #view current levels
# Relabel the levels by first creating a vector of the current levels
landcover_levels <- levels(lw_data$landcover2015_100m)
# Modify specific levels by direct replacement
landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10
# Assign the modified levels back to the factor
levels(lw_data$landcover2015_100m) <- landcover_levels
# Verify the new levels
levels(lw_data$landcover2015_100m)

#set categories of firedec
str(lw_data$ccfiredec)
lw_data$ccfiredec <- factor(lw_data$ccfiredec)
levels(lw_data$ccfiredec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(lw_data), 0.8*nrow(lw_data))
lw_train = lw_data[n,]
lw_test = lw_data[-n,] 

#########################################################################
#Run univariate model selection on environmental covariates
#Landcover
lc30_1 <- glm(presence~cclandcover2015, family=binomial, data=lw_train)
lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=lw_train)
lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=lw_train)

summary(lc100_1) #LC 100 is the best based on AIC model selection. Low, negative z scores.
#AIC - Land cover
lc_glms <- list(lc30_1, lc100_1, lc250_1)
lc_glms_names <- c('lc30_1', 'lc100_1', 'lc250_1')
aictab(cand.set = lc_glms, modnames = lc_glms_names)

#Check waterd quadratic term
lw_waterd1 <- glm(presence~waterd, family=binomial, data=lw_train)
lw_waterd2 <- glm(presence~waterd + I(waterd^2), family=binomial, data=lw_train)
# very similar delta AIC but linear term alone performed the best.  
lw_waterd_glms <- list(lw_waterd1, lw_waterd2)
lw_waterd_names <- c('lw_waterd1', 'lw_waterd2')
aictab(cand.set = lw_waterd_glms, modnames = lw_waterd_names)

#Lichen
lichbi30 <- glm(presence~cclichen_bi, family=binomial, data=lw_train)
lichenbi100 <- glm(presence~cclichen_bi100, family=binomial, data=lw_train)
lichenbi250 <- glm(presence~cclichen_bi250, family=binomial, data=lw_train)
summary(lich30_1) # Lichen bi 30 m is the best based on AIC model selection. High z scores. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~ccconifer_bi, family=binomial, data=lw_train)
conbi100 <- glm(presence~ccconifer_bi100, family=binomial, data=lw_train)
conbi250 <- glm(presence~ccconifer_bi250, family=binomial, data=lw_train)

summary(conbi250) #conbi30 is best based on AIC scores. High, positive z score. Exclude because of Landcover 250 m. 
#AIC - conifer PFT
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~ccdecid_bi, family=binomial, data=lw_train)
decidbi100 <- glm(presence~ccdecid_bi100, family=binomial, data=lw_train)
decidbi250 <- glm(presence~ccdecid_bi250, family=binomial, data=lw_train)

summary(decidbi250) #deciduous shrub bivariate 30m is the best based on AIC model selection. High, positive z score. 
#AIC - deciduous shrub PFT
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#graminoid
grambi30 <- glm(presence~ccgram_bi, family=binomial, data=lw_train)
grambi100 <- glm(presence~ccgram_bi100, family=binomial, data=lw_train)
grambi250 <- glm(presence~ccgram_bi250, family=binomial, data=lw_train)

summary(gram250_1) #graminoid binomial 30 m is the best based on AIC 
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glm_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glm_names)

#Fire
fire_1 <- glm(presence ~ ccfire2018, family = binomial, data=lw_train)
firedec_1 <- glm(presence ~ ccfiredec, family = binomial, data=lw_train)
firesum_1 <- glm(presence ~ firesum, family = binomial, data=lw_train)

summary(fire_1) #fire and firedec scored the same, fire sum was a close second. All have very low positive z scores and are not significant. 
summary(firedec_1) #will try fire and firedec in models to check if effect changes with other covariates, but may end up excluding.
summary(firesum_1) 
fire_glms <- list(fire_1, firedec_1, firesum_1)
fire_glms_names <- c('fire_1', 'firedec_1', 'firesum_1')
aictab(cand.set = fire_glms, modnames = fire_glms_names)

#Include: dem (q), slope (q), LC100, Lichenbi30, conbi30, grambi30, fire

############################################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Late Winter:
lw_data$ccslope2 <- lw_data$ccslope^2 #quadratic transformation
lw_data$ccdem2 <- lw_data$ccdem^2

#missing_cols <- setdiff(c("ccdem", "ccslope", "ccrough", "ndvi_lw", "waterd", "landcover2015_250m", "cclichen2015", "ccconifer_250", "ccdecid_bi100", "ccgram2015_250m", "cceverg2015", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018d", "distden2018_250", "lineden2018_250", "mine2018_30", "mine2018d", "road2018_30", "road2018d"), colnames(lw_data))
#missing_cols # check for mis matched names if correlation code doesn't work.


lw_dat.cor = lw_data[,c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_100m",  "cclichen_bi", "ccconifer_bi", "ccgram_bi",  "ccfire2018", "ccfiredec", "firesum", "dist2018_30",  "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")]
names(lw_dat.cor) = c("ccdem", "ccdem2", "ccslope", "ccslope2", "landcover2015_100m", "cclichen_bi", "ccconifer_bi", "ccgram_bi", "ccfire2018", "ccfiredec", "firesum", "dist2018_30", "dist2018_250", "dist2018_500", "dist2018_1000", "dist2018_2000", "dist2018_3000", "distden2018_250", "distden2018_500", "distden2018_1000", "distden2018_2000", "distden2018_3000", "lineden2018_250", "lineden2018_500", "lineden2018_1000", "lineden2018_2000", "lineden2018_3000", "mine2018_30", "mine2018_250", "mine2018_500", "mine2018_1000", "mine2018_2000", "mine2018_3000", "road2018_30", "road2018_250", "road2018_500", "road2018_1000")
cor_matrix <- cor(lw_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(lw_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/cc_distribution/output/stats/lw_correlation_matrix_core.csv", row.names = TRUE)
#slope and roughness highly correlated (0.9)
#fire variables highly correlated with each other (0.9)
#buffered road and disturbance (0.8)
#distance to mine, disturbance and roads highly correlated (0.9)
#line density and disturbance density correlated (0.8)

############################################################################
#Late Winter GLMs
#Some experimentation
#exclude distance to water
#exclude decid shrub
#check fire 2018 coefficients
#Check coefficients of other variables, keep only those that are biologically relevant 
#core data: dem (q), slope (q), LC100, Lichenbi30, conbi30, grambi30, fire
#note conifer cover and elevation are negatively correlated & need to be tested separately
#elevation mods
#all variables
lw_test1 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#no graminoid
lw_test2 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#elevation and landcover
lw_test3 <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + landcover2015_100m, family=binomial,data=lw_train)
#conifer mods
#all variables
lw_test4 <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train)

lw_test <- list(lw_test1, lw_test2, lw_test3, lw_test4, lw_test5, lw_test6)
lw_test_names <- c('lw_test1', 'lw_test2', 'lw_test3', 'lw_test4', 'lw_test5', 'lw_test6')
aictab(cand.set = lw_test, modnames = lw_test_names)

summary(lw_test4)

# elevation is the primary driver of the model
#grasslands class is not significant and has a positive effect, graminoid cover is not significant and has a negative effect
#fire has a negative effect and is close to significant. 
#remaining coefficients make sense in full model
#in conifer model, conifer cover has a significant positive effect. Makes sense biologically and agrees with landcover classes (conifer forest is reference class)
#in conifer model both grassland class and graminoid have positive, non-significant effects. 
#graminoid cover has some correlation with conifer and elevation (-0.5, 0.4). 

#Environmental base model selection
#Unique model for each season
# Phase 1: Find the best environmental base model
#Include: dem (q), slope (q), LC100, Lichenbi30, conbi30, grambi30, fire
#elevation and conifer cover are negatively correlated (-0.64)

#All top performing predictors from evaluations and fire by decade
#Elevation mods
#global
lw_glm_1e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train) #All top performing predictors from evaluations and fire sum
#no fire
lw_glm_2e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccgram_bi), family=binomial,data=lw_train) #All top performing predictors from evaluations, no fire
#no graminoid
lw_glm_3e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train) #All top performing variables from evaluations, no graminoid
#no lichen
lw_glm_4e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train) #All top performing variables from evaluations, no lichen                 
#o landcover
lw_glm_5e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train) #All top performing variables from evaluations, no landcover
#no slope
lw_glm_6e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train) #All top performing variables from evaluations, no slope
#no elevation and no conifer
lw_glm_7e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train) #all top performing variables from evaluations, no elevation
#conifer mods
#global
lw_glm_8e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#no fire
lw_glm_9e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccgram_bi), family=binomial,data=lw_train)
#no gram
lw_glm_10e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#no lichen
lw_glm_11e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(ccconifer_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#no landcover
lw_glm_12e <- glm(presence~scale(ccslope) + scale(I(ccslope^2)) + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#no slope
lw_glm_13e <- glm(presence~landcover2015_100m + scale(cclichen_bi) + scale(ccconifer_bi) + scale(ccgram_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#null model
lw_glm_14e <- glm(presence~1,family=binomial,data=lw_train) #null model

lw_e_glms <- list(lw_glm_1e, lw_glm_2e, lw_glm_3e, lw_glm_4e, lw_glm_5e, lw_glm_6e, lw_glm_7e, lw_glm_8e, lw_glm_9e, lw_glm_10e, lw_glm_11e, lw_glm_12e, lw_glm_13e, lw_glm_14e)
lw_e_glms_names <- c('lw_glm_1e', 'lw_glm_2e', 'lw_glm_3e', 'lw_glm_4e', 'lw_glm_5e', 'lw_glm_6e', 'lw_glm_7e', 'lw_glm_8e', 'lw_glm_9e', 'lw_glm_10e', 'lw_glm_11e', 'lw_glm_12e', 'lw_glm_13e', 'lw_glm_14e')
aictab(cand.set = lw_e_glms, modnames = lw_e_glms_names)
#1 e and 2e are best models; 1e has firesum (delta AIC 1.13), 2e has no fire.
#graminoid and slope are least important variables
#elevation, conifer, landcover, lichen are all important variables. 

summary(lw_glm_3e)
summary(lw_glm_11e)
summary(lw_glm_1e)
summary()

report(lw_glm_1e)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Step 1: Fit the GLM model on the training data

# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(lw_glm_3e, newdata = lw_test, type = "response")

# Step 3: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = lw_test$presence,       # Binary response variable in the test data
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
predicted_probs_test <- predict(lw_glm_3e, newdata = lw_test, type = "response")

# Step 2: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = lw_test$presence,       # Binary response variable
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

#Good calibration. model performs reasonably well. Predicted probabilities generally align with the observed rates across the 10 bins. 
#the model slightly under-predicts in the lower-middle range and slightly over-predicts for higher probability bins

#Evaluation with AUC / ROC
#unscaled model:
lw_glm_3e_us <- glm(presence~ccdem + I(ccdem^2) + ccslope + I(ccslope^2) + landcover2015_100m + cclichen_bi + ccfire2018, family=binomial,data=lw_train)
#evaluate
lw_glm_3e_eval <- pa_evaluate(predict(lw_glm_3e_us, lw_test[lw_test$presence==1, ]), predict(lw_glm_3e_us, lw_test[lw_test$presence==0, ]))
print(lw_glm_3e_eval)
plot(lw_glm_3e_eval, "ROC")
#AUC=0.849 #high AUC score, good discrimination
###########################################
#Model performance with easystats
r2(lw_glm_3e_us) #Tjur's R2 = 0.230

#Predicted occurrence map
#load rasters, resample 250 and 100 m rasters, and merge. 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccfirerast2018_core.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m_core.tif') %>%
  resample(ccmask, method="near")
#rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m2.tif') %>%
 # resample(ccmask, method="near")

landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_100m <- rast100m$landcover2015_100m |> #resample landcover 250m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(20, 10) |>
  subst(18, 16) 
cclichen_bi <- envrast$cclichen_bi
#ccconifer_bi250 <- rast250m$ccconifer_bi250 |>
  #resample(ccmask, method='near')
ccgram_bi <- envrast$ccgram_bi
ccfire2018 <- firerast$ccfire2018

rasters <- c(ccdem, ccslope, landcover2015_100m, cclichen_bi, ccgram_bi, ccfire2018) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, lw_glm_3e_us, type="response")
plot(p1) #plots!

#Add gps points to map
#filtered used points w/o covariates extracted to them:
lw_points <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Caribou Data/GPS_collar_data/cc_gps_2018_latewinter_core.csv')
st_crs(lw_points)
# convert from lat/long to UTM to match rasters
lw_points <- st_as_sf(lw_points, coords = c("location_long", "location_lat"), crs = 4326)
lw_points <- st_transform(lw_points, crs = st_crs(rasters)) #ensure crs match - doesn't seem to be converting to UTM
#check extents:
st_bbox(rasters)
st_bbox(lw_points)
plot(st_geometry(lw_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

###################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation
#add quadratic term for distance to variables.
#remember buffered road and disturbance; distance to covariates; line density and disturbance density are highly correlated and can't go in the same model. 
#model le + disturbance 30m buff
lw_glm_1d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(dist2018_30), family=binomial,data=lw_train)
#model le + disturbance 250m buff
lw_glm_2d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(dist2018_250), family=binomial,data=lw_train)
#model le + disturbance 500m buff
lw_glm_3d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(dist2018_500), family=binomial,data=lw_train)
#model 1e + disturbance 1000m buff
lw_glm_4d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(dist2018_1000), family=binomial,data=lw_train)
#dist 1000m buff * elevation
lw_glm_4di <- glm(presence~scale(ccdem) * scale(dist2018_1000) + scale(I(ccdem^2)) * scale(dist2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#model 1e + disturbance 2000m buff
lw_glm_5d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(dist2018_2000), family=binomial,data=lw_train)
#dist 2000m buff * elevation
lw_glm_5di <- glm(presence~scale(ccdem) * scale(dist2018_2000) + scale(I(ccdem^2)) * scale(dist2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#disturbance 3000 m buff
lw_glm_6d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(dist2018_3000), family=binomial, data=lw_train)
#dist 3000m buff * elevation
lw_glm_6di <- glm(presence~scale(ccdem) * scale(dist2018_3000) + scale(I(ccdem^2)) * scale(dist2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#model 1e + mining 30 m buff
lw_glm_7d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(mine2018_30), family=binomial,data=lw_train)
#model 1e + mining 250 m buff
lw_glm_8d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(mine2018_250), family=binomial,data=lw_train)
#model 1e + mining 500 m buff
lw_glm_9d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(mine2018_500), family=binomial,data=lw_train)
#model 1e + mining 1000 m buff
lw_glm_10d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(mine2018_1000), family=binomial,data=lw_train)
#mining 2000 m buff
lw_glm_11d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(mine2018_2000), family=binomial,data=lw_train)
#mining 3000 m buff
lw_glm_12d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(mine2018_3000), family=binomial,data=lw_train)
#model 1e + roads and trails 30 m buff
lw_glm_13d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(road2018_30), family=binomial,data=lw_train)
#model 1e + roads and trails 250 m buff
lw_glm_14d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(road2018_250), family=binomial,data=lw_train)
#model 1e + roads and trails 500 m buff
lw_glm_15d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(road2018_500), family=binomial,data=lw_train)
#model 1e + roads and trails 1000 m buff
lw_glm_16d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(road2018_1000), family=binomial,data=lw_train)
#road 1000 m buff * elevation
lw_glm_16di <- glm(presence~scale(ccdem) * scale(road2018_1000) + scale(I(ccdem^2)) * scale(road2018_1000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#model 1e + disturbance density 250 m radii
lw_glm_17d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(distden2018_250), family=binomial,data=lw_train)
#model 1e + disturbance density 500 m radii
lw_glm_18d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(distden2018_500), family=binomial,data=lw_train)
#model 1e + disturbance density 1000 m radii
lw_glm_19d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(distden2018_1000), family=binomial,data=lw_train)
#model 1e + disturbance density 2000 m radii
lw_glm_20d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(distden2018_2000), family=binomial,data=lw_train)
#model 1e + disturbance density 3000 m radii
lw_glm_21d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(distden2018_3000), family=binomial,data=lw_train)
#model 1e + LF density 250 m radii
lw_glm_22d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(lineden2018_250), family=binomial,data=lw_train)
#model 1e + LF density 500 m radii
lw_glm_23d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(lineden2018_500), family=binomial,data=lw_train)
#model 1e + LF density 1000 m radii
lw_glm_24d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(lineden2018_1000), family=binomial,data=lw_train)
#model 1e + LF density 2000 m radii
lw_glm_25d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(lineden2018_2000), family=binomial,data=lw_train)
#LF density 2000 m radii * elevation
lw_glm_25di <- glm(presence~scale(ccdem) * scale(lineden2018_2000) + scale(I(ccdem^2)) * scale(lineden2018_2000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#model 1e + LF density 3000 m radii
lw_glm_26d <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018) + scale(lineden2018_3000), family=binomial,data=lw_train)
#LF density 3000 m radii * elevation
lw_glm_26di <- glm(presence~scale(ccdem) * scale(lineden2018_3000) + scale(I(ccdem^2)) * scale(lineden2018_3000) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#model 1e
lw_glm_3e <- glm(presence~scale(ccdem) + scale(I(ccdem^2)) + scale(ccslope) + scale(I(ccslope^2)) + landcover2015_100m + scale(cclichen_bi) + scale(ccfire2018), family=binomial,data=lw_train)
#null model
lw_glm_27d <- glm(presence~1,family=binomial,data=lw_train)

#AIC
lw_d_glms <- list(lw_glm_1d, lw_glm_2d, lw_glm_3d, lw_glm_4d, lw_glm_4di, lw_glm_5d, lw_glm_5di, lw_glm_6d, lw_glm_6di, lw_glm_7d, lw_glm_8d, lw_glm_9d, lw_glm_10d, lw_glm_11d, lw_glm_12d, lw_glm_13d, lw_glm_14d, lw_glm_15d, lw_glm_16d, lw_glm_16di, lw_glm_17d, lw_glm_18d, lw_glm_19d, lw_glm_20d, lw_glm_21d, lw_glm_22d, lw_glm_23d, lw_glm_24d, lw_glm_25d, lw_glm_25di, lw_glm_26d, lw_glm_26di, lw_glm_27d, lw_glm_3e)
lw_d_glm_names <- c('lw_glm_1d', 'lw_glm_2d', 'lw_glm_3d', 'lw_glm_4d', 'lw_glm_4di', 'lw_glm_5d', 'lw_glm_5di', 'lw_glm_6d', 'lw_glm_6di', 'lw_glm_7d', 'lw_glm_8d', 'lw_glm_9d', 'lw_glm_10d', 'lw_glm_11d', 'lw_glm_12d', 'lw_glm_13d', 'lw_glm_14d', 'lw_glm_15d', 'lw_glm_16d', 'lw_glm_16di', 'lw_glm_17d', 'lw_glm_18d', 'lw_glm_19d', 'lw_glm_20d', 'lw_glm_21d', 'lw_glm_22d', 'lw_glm_23d', 'lw_glm_24d', 'lw_glm_25d', 'lw_glm_25di', 'lw_glm_26d', 'lw_glm_26di', 'lw_glm_27d', 'lw_glm_3e')
aictab(cand.set = lw_d_glms, modnames = lw_d_glm_names)

#disturbance with 1000 m buffer is top disuturbance variable, followed closely by line density 3000 m radii (delta AIC of 0.77)
#Disturbance density 3000 m radii close third top variable, then dist 2000 m buff. 
#all top variables are correlated, can't do a combined model.
#interaction with elevation for correlated variables improves model fit


summary(lw_glm_26di)
summary(lw_glm_4di)
summary(lw_glm_5d) #both slope terms now non-significant, graminoid non-significant
report()

## Model Calibration
#plot predicted vs observed probabilities of occurrence:

# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(lw_glm_26di, newdata = lw_test, type = "response")

# Step 3: Create a dataframe with observed outcomes and predicted probabilities
results_test <- data.frame(
  observed = lw_test$presence,       # Binary response variable in the test data
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
# Bin the predicted probabilities into deciles (10 equal-sized bins)
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
#Good calibration. model performs well. 
#the model slightly under-predicts in the lower-middle range and slightly over-predicts for higher probabilities

#Test the model
#top model with unscaled predictors
lw_glm_26di_us <- glm(presence~ccdem * lineden2018_3000 + I(ccdem^2) * lineden2018_3000 + ccslope + I(ccslope^2) + landcover2015_100m + cclichen_bi + ccfire2018, family=binomial,data=lw_train)
#evaluate / model validation
lw_glm26di_eval <- pa_evaluate(predict(lw_glm_26di_us, lw_test[lw_test$presence==1, ]), predict(lw_glm_26di_us, lw_test[lw_test$presence==0, ]))
print(lw_glm26di_eval)
plot(lw_glm26di_eval, "ROC")
#AUC=0.857  #AUC score improved from base model.
######################################
#Model performance with easystats
r2(lw_glm_26di_us) #Tjur's R2 = 0.249

windows()
check_model(lw_glm_26di_us)
#effect size plot with see pkg from easystats
plot(effectsize(lw_glm_26di_us)) +
  ggplot2::labs(title = "Clear Creek - Late Winter 2018") +
  ggplot2::scale_y_discrete(labels = c(
    "ccdem" = "Elevation",
    "lineden2018_3000" = "Line density 3000 m",
    "I(ccdem^2)" = "Elevation^2",
    "ccslope" = "Slope",
    "I(ccslope^2)" = "Slope^2",
    "landcover2015_100m5" = "Deciduous/mixed 100 m",
    "landcover2015_100m8" = "Shrublands 100 m",
    "landcover2015_100m10" = "Grasslands 100 m",
    "landcover2015_100m14" = "Wetlands 100 m",
    "landcover2015_100m16" = "Non-vegetated 100 m",
    "cclichen_bi" = "Lichen 30 m",
    "ccfire2018" = "Burns ≤ 50 yrs",
    "ccdem:lineden2018_3000" = "Elevation:Line density 3000 m",
    "lineden2018_3000:I(ccdem^2)" = "Elevation^2:Line density 3000 m"
  )) 

plot(parameters(lw_glm_26di_us)) +
  ggplot2::labs(title = "Late Winter")

#Save models
saveRDS(lw_glm_26di, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_lw_scaled.rds")
saveRDS(lw_glm_26di_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/cc_lw_unscaled.rds")

#load rasters, resample 250 and 100 m rasters, and merge. 
ccmask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/cc_core_mask_30m.tif')
ccmask <- subst(ccmask,0,NA) 
ccmask <- trim(ccmask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccenvrasters2018_core.tif')
firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccfirerast2018_core.tif')
#distbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistbuffrast2018_1.tif')
#distdenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccdistdenrast2018_1.tif')
#minebuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccminebuffrast2018_1.tif')
rast100m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_100m_core.tif') %>%
  resample(ccmask, method="near")
#rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/ccrasters2018_250m2.tif') %>%
  #resample(ccmask, method="near")
linerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/cclinedenrast2018_core.tif')

landcover_levels[landcover_levels == "20"] <- "10"  # Replace level 20 with 10
landcover_levels[landcover_levels == "18"] <- "16"  # Replace level 18 with 16
landcover_levels[landcover_levels == "15"] <- "10"  # Replace level 15 with 10

#prepare rasters
ccdem <- as.numeric(envrast$ccdem) #dem to numeric - should have been done at raster prep script
ccslope <- envrast$ccslope
landcover2015_100m <- rast100m$landcover2015_100m |> #resample landcover 250m to 30m, assign as factor, set missing category (20) to Grassland (10). Keep original variable name (even though its 30m now)
  resample(ccmask, method='near') |>
  as.factor() |>
  subst(20, 10) |>
  subst(18, 16) 
cclichen_bi <- envrast$cclichen_bi
ccgram_bi <- envrast$ccgram_bi
ccfire2018 <- firerast$ccfire2018
lineden2018_3000 <- linerast$lineden2018_3000

rasters <- c(ccdem, ccslope, landcover2015_100m, cclichen_bi, ccgram_bi, ccfire2018, lineden2018_3000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, lw_glm_26di_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/cc_lw_mod26di_prediction.tif')

#Add lw used points (all) to map:
plot(st_geometry(lw_points), add = TRUE, col = "red", pch = 19, cex = 0.6) 

####################################################
####Calculate peak of quadratic relationship for elevation
windows()
summary(lw_glm_26di_us) #get unscaled coefficients for elevation
#dem = 0.06837, dem2 = -3.119*10^-5
lw_dem_peak = -(0.06837 / (2 * (-3.119*10^-5)))
lw_dem_peak # peak occurs at 1096.024 metres 

#Plot Relationship
# Coefficients from the model
beta_0 <- -39.19          # Intercept
beta_1 <- 0.06837         # Coefficient for dem
beta_2 <- -3.119e-05      # Coefficient for I(dem^2)

# Create a sequence of dem values around the observed range
dem_values <- seq(500, 2000, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * dem_values + beta_2 * dem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(dem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Late Winter")
abline(v = 1096.024, col = "red", lty = 2)  # Add vertical line at the peak
