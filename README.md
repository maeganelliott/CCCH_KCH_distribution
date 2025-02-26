## CCCH_KCH_distribution
 Public respository for R code used to assess the effects of landscape change on the distributions of the Clear Creek and Klaza caribou herds in the Yukon, Canada, for M. Elliott's MSc Thesis.
 The Thesis is available here: 

#####################################################
Overview:
This repository does not contain the data used because it is proprietary. The analysis used caribou location data from satellite/Global Positioning System (GPS) caribou collars and Very High Frequency (VHF) caribou collars. The analysis was split into separate time stamps for each herd based on the collar data, and each time stamp was further divided into seasons. This repository only includes the scripts used to run the Resource Selection Functions (RSFs), which were Generalized Linear Models (GLMs). A breakdown of the analysis steps from data preparation and cleaning to model selection is provided in *Flow_chart.pdf*.

I used the renv (https://rstudio.github.io/renv/index.html) R package to take a snapshot of my R scripts and the associated versions of the packages I used to ensure they will run later if needed. Follow the package instructions; call renv::restore() to use the same package versions required by the project that are stored in the lockfile. 

#####################################################
Explanation of Scripts:
Each caribou herd has its own set of scripts stored in the ClearCreek_RSFs and Klaza_RSFs folders. 
01_gps_glm_mods_ew = GPS collar data GLMs for the early winter season
02_gps_glm_mods_lw = GPS collar data GLMs for the late winter season
03_gps_glm_mods_s = GPS collar data GLMs for the summer season
04_gps_glm_mods_f = GPS collar data GLMs for the fall-rut season
05_vhf_glm_mods = VHF collar data GLMs for the snow and snow-free seasons
