## Clear Creek and Klaza Caribou Herd Distribution Changes

 Public repository for supplementary material for M. Elliott's MSc Thesis and corresponding journal article. This includes:
 - R code used to assess the effects of landscape change on the distributions of the Clear Creek (CCCH) and Klaza (KCH) caribou herds in the Yukon, Canada
 - Community knowledge based definitions of factors affecting each caribou herd
 - A summary of the number of collared caribou and locations per season included in RSFs
 - A description of RSF covariates
 - AICc tables from the model selection process
 - RSF Regression tables

M. Elliott's MSc Thesis can be accessed here: https://ualberta.scholaris.ca/items/c2591947-7bc7-4f91-b9e1-b2f5a7f9aac8 

###########################################################################################################
### Overview:

This public repository does not contain the data used because it is proprietary. Resource Selection Functions (RSFs) used caribou location data from satellite/Global Positioning System (GPS) caribou collars and Very High Frequency (VHF) caribou collars. The analysis was split into separate time stamps for each herd based on the collar data, and each time stamp was further divided into seasons. This repository only includes the scripts used to run the Resource Selection Functions (RSFs), which were Generalized Linear Models (GLMs). A breakdown of the analysis steps from data preparation and cleaning to model selection is provided in *Flow_chart.pdf*.

#########################################################################################################
### Explanation of Scripts:

Each caribou herd has its own set of scripts stored in the ClearCreek_RSFs and Klaza_RSFs folders. 

01_gps_glm_mods_ew = GPS collar data GLMs for the early winter season

02_gps_glm_mods_lw = GPS collar data GLMs for the late winter season

03_gps_glm_mods_s = GPS collar data GLMs for the summer season

04_gps_glm_mods_f = GPS collar data GLMs for the fall-rut season

05_vhf_glm_mods = VHF collar data GLMs for the snow and snow-free seasons

#######################################################################################################
### Temporal Disturbance Database:
Please visit https://github.com/beaconsproject/disturbance_mapping for the temporal disturbance database used in this analysis (under Master's thesis), as well as a number of other disturbance mapping resources. 
