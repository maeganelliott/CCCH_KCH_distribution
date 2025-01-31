#Script for Klaza gps glm models
## LATE WINTER

#https://github.com/rygill/caribou_movement_ecology/blob/main/scripts/07.HR_Regression.r

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
library(dplyr)
library(ggplot2)
library(car)
library(AICcmodavg)

set.seed(123)

#Import used and available points with covariates previously extracted to them
lw_data <- read.csv('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_lw_data3.csv')

#Ensure categorical variables are set as factors (categories) -- categorical variables don't need to be scaled. 
str(lw_data$landcover2015_250m)
lw_data$landcover2015_250m <- factor(lw_data$landcover2015_250m)
levels(lw_data$landcover2015_250m) #view current levels

#set categories of firedec
str(lw_data$firedec)
lw_data$firedec <- factor(lw_data$firedec)
levels(lw_data$firedec)

#Split into Training & Testing datasets 80% - 20%
n = sample(nrow(lw_data), 0.8*nrow(lw_data))
lw_train = lw_data[n,]
lw_test = lw_data[-n,] 

###########################################################################
#Univariate Model Selection
#Run univariate model selection on covariates that were unclear from above evaluations
#waterd
lw_waterd1 <- glm(presence~waterd, family=binomial, data=lw_train)
lw_waterd2 <- glm(presence~waterd + I(waterd^2), family=binomial, data=lw_train) #quadratic version performs better. 

lw_waterd_glms <- list(lw_waterd1, lw_waterd2)
lw_waterd_glms_names <- c('lw_waterd1', 'lw_waterd2')
aictab(cand.set = lw_waterd_glms, modnames = lw_waterd_glms_names)

#Landcover
lc30_1 <- glm(presence~landcover2015, family=binomial, data=lw_train)
lc100_1 <- glm(presence~landcover2015_100m, family=binomial, data=lw_train)
lc250_1 <- glm(presence~landcover2015_250m, family=binomial, data=lw_train)

summary(lc250_1) #LC 250 is the best based on AIC model selection. Low, negative z scores.
#AIC - Land cover
lc_glms <- list(lc30_1, lc100_1, lc250_1)
lc_glms_names <- c('lc30_1', 'lc100_1', 'lc250_1')
aictab(cand.set = lc_glms, modnames = lc_glms_names)

#Lichen
lichbi30 <- glm(presence~lichen_bi, family=binomial, data=lw_train)
lichenbi100 <- glm(presence~lichen_bi100, family=binomial, data=lw_train)
lichenbi250 <- glm(presence~lichen_bi250, family=binomial, data=lw_train)
summary(lich30_1) # Lichen bi 30 m is the best based on AIC model selection. High z scores. 
#AIC - lichen PFT
lich_glms <- list(lichbi30, lichenbi100, lichenbi250)
lich_glms_names <- c('lichbi30', 'lichenbi100', 'lichenbi250')
aictab(cand.set = lich_glms, modnames = lich_glms_names)

#Conifer
conbi30 <- glm(presence~conifer_bi, family=binomial, data=lw_train)
conbi100 <- glm(presence~conifer_bi100, family=binomial, data=lw_train)
conbi250 <- glm(presence~conifer_bi250, family=binomial, data=lw_train)

summary(conbi250) #conbi250 is best based on AIC scores. High, positive z score. Exclude because of Landcover 250 m. 
#AIC - conifer PFT
con_glms <- list(conbi30, conbi100, conbi250)
con_glms_names <- c('conbi30', 'conbi100', 'conbi250')
aictab(cand.set = con_glms, modnames = con_glms_names)

#Deciduous shrub
decidbi30 <- glm(presence~decid_bi, family=binomial, data=lw_train)
decidbi100 <- glm(presence~decid_bi100, family=binomial, data=lw_train)
decidbi250 <- glm(presence~decid_bi250, family=binomial, data=lw_train)

summary(decidbi30) #deciduous shrub bivariate 30m is the best based on AIC model selection. High, positive z score. 
#AIC - deciduous shrub PFT
decid_glms <- list(decidbi30, decidbi100, decidbi250)
decid_glms_names <- c('decidbi30', 'decidbi100', 'decidbi250')
aictab(cand.set = decid_glms, modnames = decid_glms_names)

#graminoid
grambi30 <- glm(presence~gram_bi, family=binomial, data=lw_train)
grambi100 <- glm(presence~gram_bi100, family=binomial, data=lw_train)
grambi250 <- glm(presence~gram_bi250, family=binomial, data=lw_train)

summary(gram250_1) #graminoid binomial 250 m is the best based on AIC - exclude because of landcover 250 m
gram_glms <- list(grambi30, grambi100, grambi250)
gram_glm_names <- c('grambi30', 'grambi100', 'grambi250')
aictab(cand.set = gram_glms, modnames = gram_glm_names)


#Fire
fire_1 <- glm(presence ~ fire2015, family = binomial, data=lw_train)
firedec_1 <- glm(presence ~ firedec, family = binomial, data=lw_train)
firesum_1 <- glm(presence ~ firesum, family = binomial, data=lw_train)

summary(fire_1) #firesum performed the best. Will still test other fire variables to confirm best fit.
summary(firedec_1) 
summary(firesum_1) 
fire_glms <- list(fire_1, firedec_1, firesum_1)
fire_glms_names <- c('fire_1', 'firedec_1', 'firesum_1')
aictab(cand.set = fire_glms, modnames = fire_glms_names)

#Include: dem (q), slope (q), LC250, Lichenbi 30, coniferbi 250, grambi 250, fire
#excluding waterd because no relationship based on distribution and primary  water courses probably not relevent in winter
#excluding deciduous shrub because not biologically important in winter

############################################################################################
#Check for collinearity among predictors using corrplot package and function (from Gill repository)
#Late Winter:
lw_data$slope2 <- lw_data$slope^2 #quadratic transformation
lw_data$kdem2 <- lw_data$kdem^2

lw_dat.cor = lw_data[,c("kdem", "kdem2", "slope", "slope2", "landcover2015_250m", "lichen_bi", "conifer_bi250", "decid_bi", "gram_bi250",  "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000")]
names(lw_dat.cor) = c("kdem", "kdem2", "slope", "slope2", "landcover2015_250m", "lichen_bi", "conifer_bi250", "decid_bi", "gram_bi250", "fire2015", "firedec", "firesum", "dist2015_30", "dist2015_250", "dist2015_500", "dist2015_1000", "dist2015_2000", "dist2015_3000", "distden2015_250", "distden2015_500",  "distden2015_1000", "distden2015_2000", "distden2015_3000", "lineden2015_250", "lineden2015_500", "lineden2015_1000", "lineden2015_2000", "lineden2015_3000", "mine2015_30", "mine2015_250", "mine2015_500", "mine2015_1000", "mine2015_2000", "mine2015_3000", "road2015_30", "road2015_250", "road2015_500", "road2015_1000")
cor_matrix <- cor(lw_dat.cor, use = "pairwise.complete.obs")
print(cor_matrix)
windows()
corrplot(cor(lw_dat.cor), method = 'number') 
write.csv(cor_matrix, file = "C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/stats/lw_correlation_matrix.csv", row.names = TRUE)

##############################################################################
## LATE WINTER GLM MODEL SELECTION ##

#Initial model exploration
lw_test1 <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train)
#landcover only
lw_test2 <- glm(presence~landcover2015_250m, family=binomial,data=lw_train)

summary(lw_test2)
#In full model only the non-vegetated landcover class is not significant
#coefficients for each variable generally make sense for the season
#elevation has a larger effect than I would expect for this season
#proceed with all variables to AIC model selection


#Starting with GLM models first
#Unique model for each season
# Phase 1: Find the best environmental base model
#Include: dem (q), slope (q), LC250, Lichenbi 30, coniferbi 250, grambi 250, fire
#Starting with glm (no mixed effect), standardized predictors:
#All top performing predictors from evaluations and fire
lw_glm_1e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train) #All top performing predictors from evaluations and fire by decade
#All top performing predictors from evaluations, no fire
lw_glm_2e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250), family=binomial,data=lw_train) #All top performing predictors from evaluations, no fire
#All top performing variables from evaluations, no graminoid
lw_glm_3e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(fire2015), family=binomial,data=lw_train) #All top performing variables from evaluations, no graminoid
#All top performing variables from evaluations, no conifer
lw_glm_4e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train) #All top performing variables from evaluations, no conifer
#All top performing variables from evaluations, no lichen
lw_glm_5e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train) #All top performing variables from evaluations, no lichen                 
#All top performing variables from evaluations, no landcover
lw_glm_6e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train) #All top performing variables from evaluations, no landcover
#All top performing variables from evaluations, no slope
lw_glm_7e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train) #All top performing variables from evaluations, no slope
#all top performing variables from evaluations, no elevation
lw_glm_8e <- glm(presence~scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train) #all top performing variables from evaluations, no elevation
#null model
lw_glm_9e <- glm(presence~1,family=binomial,data=lw_train) #null model

lw_e_glms <- list(lw_glm_1e, lw_glm_2e, lw_glm_3e, lw_glm_4e, lw_glm_5e, lw_glm_6e, lw_glm_7e, lw_glm_8e, lw_glm_9e)
lw_e_glms_names <- c('lw_glm_1e', 'lw_glm_2e', 'lw_glm_3e', 'lw_glm_4e', 'lw_glm_5e', 'lw_glm_6e', 'lw_glm_7e', 'lw_glm_8e', 'lw_glm_9e')
aictab(cand.set = lw_e_glms, modnames = lw_e_glms_names)
#full model has the best fit
#graminoid and conifer are least important variables
#elevation, landcover and fire are most important variables. 
summary(lw_glm_1e)
summary()

report(lw_glm_1e)

## Model Calibration
#plot predicted vs observed probabilities of occurrence:
# Step 1: Fit the GLM model on the training data (already done in your case with lw_train)
# Assuming the model is already fit and called lw_glm_1e

# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(lw_glm_1e, newdata = lw_test, type = "response")

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
predicted_probs_test <- predict(lw_glm_1e, newdata = lw_test, type = "response")

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

#Evaluation with AUC / ROC
#unscaled model:
lw_glm_1e_us <- glm(presence~kdem + I(kdem^2) + slope + I(slope^2) + landcover2015_250m + lichen_bi + conifer_bi250 + gram_bi250 + fire2015, family=binomial,data=lw_train)
#evaluate
lw_glm_1e_eval <- pa_evaluate(predict(lw_glm_1e_us, lw_test[lw_test$presence==1, ]), predict(lw_glm_1e_us, lw_test[lw_test$presence==0, ]))
print(lw_glm_1e_eval)
plot(lw_glm_1e_eval, "ROC")
#AUC=0.881 #high AUC score, good discrimination
#################
#Model performance with easystats
r2(lw_glm_1e_us) #Tjur's R2 = 0.298


#predicted occurrence map
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kfirerast2015_1.tif')
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast$kdem) #dem to numeric - should have been done at raster prep script
slope <- envrast$slope
landcover2015_250m <- as.factor(rast250m$landcover2015_250m) 
lichen_bi <- envrast$lichen_bi
conifer_bi250 <- rast250m$conifer_bi250 |>
  resample(mask, method='near')
gram_bi250 <- rast250m$gram_bi250 |>
  resample(mask, method='near')
fire2015 <- firerast$fire2015

rasters <- c(kdem, slope, landcover2015_250m, lichen_bi, conifer_bi250, gram_bi250, fire2015) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, lw_glm_1e_us, type="response")
plot(p1) #plots!

#Add gps points to map
#filtered used points w/o covariates extracted to them:
lw_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_lw_daily.shp')

#Plot the raster and add points on top
plot(st_geometry(lw_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

###################################################################################
# Phase 2: Test effect of different disturbance variables on base environmental model. 
#only test the variables that make sense based on previous evaluation
#add quadratic term for distance to variables.
#remember buffered road and disturbance; distance to covariates; line density and disturbance density are highly correlated and can't go in the same model. 
#model le + disturbance 30m buff
lw_glm_1d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(dist2015_30), family=binomial,data=lw_train)
#model le + disturbance 250m buff
lw_glm_2d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(dist2015_250), family=binomial,data=lw_train)
#model le + disturbance 500m buff
lw_glm_3d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(dist2015_500), family=binomial,data=lw_train)
#model 1e + disturbance 1000m buff
lw_glm_4d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(dist2015_1000), family=binomial,data=lw_train)
#model 1e + disturbance 2000m buff
lw_glm_5d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(dist2015_2000), family=binomial,data=lw_train)
#model 1e + disturbance 3000 m buff
lw_glm_6d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(dist2015_3000), family=binomial,data=lw_train)
#model 1e + mining 30 m buff
lw_glm_7d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(mine2015_30), family=binomial,data=lw_train)
#model 1e + mining 250 m buff
lw_glm_8d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(mine2015_250), family=binomial,data=lw_train)
#model 1e + mining 500 m buff
lw_glm_9d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(mine2015_500), family=binomial,data=lw_train)
#model 1e + mining 1000 m buff
lw_glm_10d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(mine2015_1000), family=binomial,data=lw_train)
#model 1e + mining 2000 m buff
lw_glm_11d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(mine2015_2000), family=binomial,data=lw_train)
#model 1e + mining 3000 m buff
lw_glm_12d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(mine2015_3000), family=binomial,data=lw_train)
#model 1e + roads and trails 30 m buff
lw_glm_13d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(road2015_30), family=binomial,data=lw_train)
#model 1e + roads and trails 250 m buff
lw_glm_14d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(road2015_250), family=binomial,data=lw_train)
#model 1e + roads and trails 500 m buff
lw_glm_15d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(road2015_500), family=binomial,data=lw_train)
#model 1e + roads and trails 1000 m buff
lw_glm_16d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(road2015_1000), family=binomial,data=lw_train)
#model 1e + disturbance density 250 m radii
lw_glm_17d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(distden2015_250), family=binomial,data=lw_train)
#model 1e + disturbance density 500 m radii
lw_glm_18d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(distden2015_500), family=binomial,data=lw_train)
#model 1e + disturbance density 1000 m radii
lw_glm_19d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(distden2015_1000), family=binomial,data=lw_train)
#model 1e + disturbance density 2000 m radii
lw_glm_20d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(distden2015_2000), family=binomial,data=lw_train)
#model 1e + disturbance density 3000 m radii
lw_glm_21d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(distden2015_3000), family=binomial,data=lw_train)
#model 1e + LF density 250 m radii
lw_glm_22d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(lineden2015_250), family=binomial,data=lw_train)
#model 1e + LF density 500 m radii
lw_glm_23d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(lineden2015_500), family=binomial,data=lw_train)
#model 1e + LF density 1000 m radii
lw_glm_24d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(lineden2015_1000), family=binomial,data=lw_train)
#model 1e + LF density 2000 m radii
lw_glm_25d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(lineden2015_2000), family=binomial,data=lw_train)
#model 1e + LF density 3000 m radii
lw_glm_26d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(lineden2015_3000), family=binomial,data=lw_train)
#model 1e
lw_glm_1e <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015), family=binomial,data=lw_train)
#null model
lw_glm_27d <- glm(presence~1,family=binomial,data=lw_train)
#model 1e + lineden 3000 + dist 3000 m buff
lw_glm_28d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(lineden2015_3000) + scale(dist2015_3000), family=binomial,data=lw_train)
#model 1e + distden 3000 + dist 3000 m buff
lw_glm_29d <- glm(presence~scale(kdem) + scale(I(kdem^2)) + scale(slope) + scale(I(slope^2)) + landcover2015_250m + scale(lichen_bi) + scale(conifer_bi250) + scale(gram_bi250) + scale(fire2015) + scale(distden2015_3000) + scale(dist2015_3000), family=binomial,data=lw_train)

#AIC
lw_d_glms <- list(lw_glm_1d, lw_glm_2d, lw_glm_3d, lw_glm_4d, lw_glm_5d, lw_glm_6d, lw_glm_7d, lw_glm_8d, lw_glm_9d, lw_glm_10d, lw_glm_11d, lw_glm_12d, lw_glm_13d, lw_glm_14d, lw_glm_15d, lw_glm_16d, lw_glm_17d, lw_glm_18d, lw_glm_19d, lw_glm_20d, lw_glm_21d, lw_glm_22d, lw_glm_23d, lw_glm_24d, lw_glm_25d, lw_glm_26d, lw_glm_27d, lw_glm_28d, lw_glm_29d, lw_glm_1e)
lw_d_glm_names <- c('lw_glm_1d', 'lw_glm_2d', 'lw_glm_3d', 'lw_glm_4d', 'lw_glm_5d', 'lw_glm_6d', 'lw_glm_7d', 'lw_glm_8d', 'lw_glm_9d', 'lw_glm_10d', 'lw_glm_11d', 'lw_glm_12d', 'lw_glm_13d', 'lw_glm_14d', 'lw_glm_15d', 'lw_glm_16d', 'lw_glm_17d', 'lw_glm_18d', 'lw_glm_19d', 'lw_glm_20d', 'lw_glm_21d', 'lw_glm_22d', 'lw_glm_23d', 'lw_glm_24d', 'lw_glm_25d', 'lw_glm_26d', 'lw_glm_27d', 'lw_glm_28d', 'lw_glm_29d', 'lw_glm_1e')
aictab(cand.set = lw_d_glms, modnames = lw_d_glm_names)

#larger density radii and larger buffer sizes performed the best
# linedenity 3000 top variable, followed by distden 3000
#dist 3000 m buff top buffered variable
#most disturbance variables have significant negative effects

#lineden and distden 3000 are correlated
#lineden and dist 3000 m buff can go in same model
#distden 3000 and dist 3000 m buff can go in same model


summary(lw_glm_6d)
summary(lw_glm_28d)

report()

## Model Calibration
#plot predicted vs observed probabilities of occurrence:

# Step 2: Predict probabilities on the testing data
predicted_probs_test <- predict(lw_glm_28d, newdata = lw_test, type = "response")

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
#Good calibration. model performs very well. Deviation increased slightly >0.7 predicted probabilities.
#the model slightly under-predicts or over predicts across the range of predicted probabilities 


#Test the model
#top model with unscaled predictors
lw_glm_28d_us <- glm(presence~kdem + I(kdem^2) + slope + I(slope^2) + landcover2015_250m + lichen_bi + conifer_bi250 + gram_bi250 + fire2015 + lineden2015_3000 + dist2015_3000, family=binomial,data=lw_train)

#evaluate / model validation
lw_glm28d_eval <- pa_evaluate(predict(lw_glm_28d_us, lw_test[lw_test$presence==1, ]), predict(lw_glm_28d_us, lw_test[lw_test$presence==0, ]))
print(lw_glm28d_eval)
plot(lw_glm28d_eval, "ROC")
#AUC=0.901  #excellent discrimination, AUC score improved from base model.
#################
#Model performance with easystats
r2(lw_glm_28d_us) #Tjur's R2 = 0.339

windows()
check_model(lw_glm_28d_us)
#effect size plot with see pkg from easystats
plot(effectsize(lw_glm_28d_us)) +
  ggplot2::labs(title = "Klaza - Late Winter 2015") +
  ggplot2::scale_y_discrete(labels = c(
    "kdem" = "Elevation",
    "I(kdem^2)" = "Elevation^2",
    "slope" = "Slope",
    "I(slope^2)" = "Slope^2",
    "landcover2015_250m5" = "Deciduous/mixed 250 m",
    "landcover2015_250m8" = "Shrublands 250 m",
    "landcover2015_250m10" = "Grasslands 250 m",
    "landcover2015_250m14" = "Wetlands 250 m",
    "landcover2015_250m16" = "Non-vegetated 250 m",
    "lichen_bi" = "Lichen 30 m",
    "conifer_bi250" = "Conifer 250 m",
    "gram_bi250" = "Graminoid 250 m",
    "fire2015" = "Burns ≤ 50 yrs",
    "lineden2015_3000" = "Line density 3000 m",
    "dist2015_3000" = "Disturbance 3000 m"
  )) 

plot(parameters(lw_glm_28d_us)) +
  ggplot2::labs(title = "Late Winter")

#Save models
saveRDS(lw_glm_28d, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_lw_scaled.rds")
saveRDS(lw_glm_28d_us, file = "G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/GLM_mods/k_lw_unscaled.rds")

#load rasters, resample 250 and 100 m rasters, and merge. 
mask <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/k_mask_5k_30m.tif')
mask <- subst(mask,0,NA) 
mask <- trim(mask)
envrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kenvrasters2015_2.tif')
firerast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kfirerast2015_1.tif')
distbuffrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/kdistbuffrast2015_1.tif')
linedenrast <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/klinedenrast2015_1.tif')
rast250m <- rast('G:/Shared drives/Clear Creek and Klaza Data/Predictor Rasters/krasters2015_250m2.tif') %>%
  resample(mask, method="near")

#prepare rasters
kdem <- as.numeric(envrast$kdem) #dem to numeric - should have been done at raster prep script
slope <- envrast$slope
landcover2015_250m <- as.factor(rast250m$landcover2015_250m) 
lichen_bi <- envrast$lichen_bi
conifer_bi250 <- rast250m$conifer_bi250 |>
  resample(mask, method='near')
gram_bi250 <- rast250m$gram_bi250 |>
  resample(mask, method='near')
fire2015 <- firerast$fire2015
lineden2015_3000 <- linedenrast$lineden2015_3000
dist2015_3000 <- distbuffrast$dist2015_3000

rasters <- c(kdem, slope, landcover2015_250m, lichen_bi, conifer_bi250, gram_bi250, fire2015, lineden2015_3000, dist2015_3000) #combine rasters. Make sure they are separated by commas.
#run predict function
p1 <- predict(rasters, lw_glm_28d_us, type="response")
plot(p1) #plots!
writeRaster(p1, 'G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_lw_mod28d_prediction_3.tif')

#reload prediction raster
#lw_prediction <- rast('G:/Shared drives/Clear Creek and Klaza Data/Analysis Telem/glm_rasters/k_lw_mod28d_prediction.tif')
#Add gps points to map
#filtered used points w/o covariates extracted to them:
lw_points <- st_read('C:/Users/info/Documents/Github_respositories/NMC_distribution/k_distribution/output/data_cleaned/kgps_lw_daily.shp')

#Plot the raster and add points on top
plot(st_geometry(lw_points), add = TRUE, col = "red", pch = 19, cex = 0.6)

###Calculate peak of quadratic relationship for elevation

summary(lw_glm_28d_us) #get unscaled coefficients for kdem & kdem2
lw_dem_peak = -(0.06493 / (2 * (-2.839*10^-5)))
lw_dem_peak # peak occurs at 1143.536 metres 

#Plot Relationship
# Coefficients from the model
beta_0 <- -39.65          # Intercept
beta_1 <- 0.06493        # Coefficient for kdem
beta_2 <- -2.839e-05      # Coefficient for I(kdem^2)

# Create a sequence of kdem values around the observed range
kdem_values <- seq(500, 2500, by = 10)  # Adjust range based on dataset

# Calculate logit values
logit_values <- beta_0 + beta_1 * kdem_values + beta_2 * kdem_values^2

# Convert logit to probability
prob_values <- 1 / (1 + exp(-logit_values))

# Plot the relationship
plot(kdem_values, prob_values, type = "l", col = "blue", lwd = 2,
     xlab = "kdem", ylab = "Probability of Presence",
     main = "Quadratic Relationship between Elevation and Presence in Late Winter")
abline(v = 1143.54, col = "red", lty = 2)  # Add vertical line at the peak
