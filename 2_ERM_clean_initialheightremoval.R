###########################################################
####### PREPARATION OF RAW DATA, CALCULATION OF MGR #######
###########################################################

####### PURPOSE OF THIS SCRIPT:
# predict aboveground biomass accounting for initial height (via model prediction)

####### RAW DATA USED:
# ERM_data_abovegroundind.rds
# ERM_data_abovegroundmeso.rds

####### DATA CREATED:
# ERM_data_abovegroundind_predicted.rds --> contains *predicted* aboveground individual plant cleaned data

###### PLOTS CREATED:
# NA

#######################
#### Load Packages ####
#######################

library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)

#######################
###### Read Data ######
#######################

dat <- readRDS("data/ERM_data_abovegroundind.rds") %>% rename(AGB = above_biomass_g, hti = height_i_cm) %>%
  drop_na("AGB", "hti", "p_level","f_level","block", "meso_num")

#### Ind biomass model ####
m1 <- lmer(log(AGB) ~ p_level*f_level + hti + (1|block), dat = dat)
anova(m1)

#### Generate Ind predictions ####

dat.pred <- dat %>% mutate(AGB_pred = exp(predict(m1)))
saveRDS(dat.pred, "data/ERM_data_abovegroundind_predicted.rds")

#### Group at mesocosm level ####
####### Join with metadata ######
meso.a <- readRDS("data/ERM_data_abovegroundmeso.rds") %>% select(-aboveBEF_biomass_g)

meso.a.pred <- dat.pred %>%
  group_by(meso_num) %>%
  summarise(aboveBEF_biomass_g = sum(AGB_pred)) %>%
  left_join(meso.a, by = "meso_num")

#######################
######## Save #########
#######################

saveRDS(meso.a.pred, "data/ERM_data_abovegroundmeso_predicted.rds")
