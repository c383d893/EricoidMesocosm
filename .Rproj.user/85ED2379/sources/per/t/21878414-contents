###########################################################
####### PREPARATION OF RAW DATA, CALCULATION OF MGR #######
###########################################################

####### PURPOSE OF THIS SCRIPT:
# clean raw biomass data and calculate MGR
# for fungal species and fungal level

####### RAW DATA USED:
# ERM_data_abovegroundind.rds

####### DATA CREATED:
# MGR_biomass_incl_fspecies_prep.rds --> contains cleaned data with a new column containing MGR
# MGR_biomass_prep.rds --> contains cleaned data with a new column containing MGR

###### PLOTS CREATED:
# NA

#######################
#### Load Packages ####
#######################

library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(multcomp)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(svglite)

############################
###### Fungal species ######
############################

#######################
###### Read Data ######
#######################

### load data
ERM_data <- readRDS("data/ERM_data_abovegroundind.rds")

### only keep rows where f_level == 0 or 1
ERM_data <- ERM_data %>%
  filter(f_level %in% c(0, 1))

### remove all rows where there is an NA in above_biomass_g
ERM_data <- ERM_data[!is.na(ERM_data$above_biomass_g),]

### remove unnecessary columns 
ERM_data <- ERM_data %>% dplyr::select(-c("belowBEF_biomass_g", "height_f_cm", "dead"))

#############################################
#### response based on abvg dry biomass  ####
#############################################

###### group data per plant & calculate mean final biomass for all control samples (sterile samples) per plant species

# sort data by plant species
ERM_data <- ERM_data[order(ERM_data$plant_sp), ]

# convert all columns except for biomass to a factor (so it's not a continuous variable where the mean is calculated afterwards)
ERM_data$meso_num <- as.factor(ERM_data$meso_num)
ERM_data$plant_id <- as.factor(ERM_data$plant_id)
ERM_data$plant_sp <- as.factor(ERM_data$plant_sp)
ERM_data$block <- as.factor(ERM_data$block)
ERM_data$p_level <- as.factor(ERM_data$p_level)
ERM_data$f_level <- as.factor(ERM_data$f_level)

###### select only rows of the control and calculate mean per plant species per diversity

# Ensure abvg_dry_biomass is numeric
ERM_data <- ERM_data %>%
  mutate(above_biomass_g = as.numeric(above_biomass_g)) 

# Calculate mean for f_level == 0 per plant_sp and p_level
mean_df <- ERM_data %>%
  filter(f_level == 0) %>% 
  group_by(plant_sp, p_level) %>% 
  summarise(mean_str = round(mean(above_biomass_g, na.rm = TRUE), 10), .groups = "drop") 

# Merge the calculated means back to the original dataframe if needed
ERM_data <- ERM_data %>%
  left_join(mean_df, by = c("plant_sp", "p_level"))

# exclude PJA p_level = 2 (because no data on PJA in f_level = 0 and p_level = 2 --> not possible to calculate MGR)
ERM_data <- ERM_data %>%
  filter(mean_str != 0)

####### calculate response of each sample normalized to mean of the sterile samples for the same plant species at the same diversity level
ERM_data <- ERM_data %>%
  mutate(
    response = if_else(
      f_level != 0 & !is.na(mean_str), 
      (above_biomass_g - mean_str) / mean_str, 
      NA_real_ # Assign NA for cases where f_level == 0 or mean_str is missing
    )
  )

#######################
######## Save #########
#######################

saveRDS(ERM_data, file = "data/MGR_biomass_incl_fspecies_prep.rds")

############################
####### Fungal level #######
############################

#######################
###### Read Data ######
#######################

###### load data
ERM_data <- readRDS("data/ERM_data_abovegroundind.rds")

###### remove all rows where there is an NA in above_biomass_g
ERM_data <- ERM_data[!is.na(ERM_data$above_biomass_g),]

###### remove unnecessary columns 
ERM_data <- ERM_data %>% dplyr::select(-c("belowBEF_biomass_g", "height_f_cm", "dead"))

#############################################
#### response based on abvg dry biomass  ####
#############################################

###### group data per plant & calculate mean final biomass for all control samples (sterile samples) per plant species

# sort data by plant species
ERM_data <- ERM_data[order(ERM_data$plant_sp), ]

# convert all columns except for biomass to a factor 
ERM_data$meso_num <- as.factor(ERM_data$meso_num)
ERM_data$plant_id <- as.factor(ERM_data$plant_id)
ERM_data$plant_sp <- as.factor(ERM_data$plant_sp)
ERM_data$block <- as.factor(ERM_data$block)
ERM_data$p_level <- as.factor(ERM_data$p_level)
ERM_data$f_level <- as.factor(ERM_data$f_level)

###### select only rows of the control and calculate mean per plant species per diversity

# make sure abvg_dry_biomass is numeric
ERM_data <- ERM_data %>%
  mutate(above_biomass_g = as.numeric(above_biomass_g)) 

# calculate mean for f_level == 0 per plant_sp and p_level
mean_df <- ERM_data %>%
  filter(f_level == 0) %>% 
  group_by(plant_sp, p_level) %>% 
  summarise(mean_str = round(mean(above_biomass_g, na.rm = TRUE), 10), .groups = "drop") 

# merge the calculated means back to the original dataframe if needed
ERM_data <- ERM_data %>%
  left_join(mean_df, by = c("plant_sp", "p_level"))

# exclude PJA p_level = 2 (because no data on PJA in f_level = 0 and p_level = 2 --> not possible to calculate MGR)
ERM_data <- ERM_data %>%
  filter(mean_str != 0)

####### calculate response of each sample normalized to mean of the sterile samples for the same plant species at the same diversity level
ERM_data <- ERM_data %>%
  mutate(
    response = if_else(
      f_level != 0 & !is.na(mean_str), 
      (above_biomass_g - mean_str) / mean_str, 
      NA_real_ # Assign NA for cases where f_level == 0 or mean_str is missing
    )
  )

#######################
######## Save #########
#######################

saveRDS(ERM_data, file = "data/MGR_biomass_prep.rds")
