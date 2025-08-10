###########################################################
############# MAIN CLEANING OF MESOCOSM DATA ##############
###########################################################

####### PURPOSE OF THIS SCRIPT:
# clean raw biomass data
# visualize data

####### RAW DATA USED:
# Ericoid_mesocosm_combinations.csv
# ERM_dataclean_April42025.txt

####### DATA CREATED:
# MGR_biomass_prep.csv --> contains cleaned data with a new column containing MGR
#ERM_data_abovegroundmeso.rds  --> contains aboveground mesocosm cleaned data
#ERM_data_belowgroundmeso.rds  --> contains belowground mesocosm cleaned data
#ERM_data_abovegroundind.rds  --> contains aboveground individual plant cleaned data

###### PLOTS CREATED:
# ERM_aboveground_ind_data.jpg
# ERM_aboveground_mesocosm_data.jpg
# ERM_belowground_mesocosm_data.jpg

#######################
#### Load Packages ####
#######################

library(tidyverse)
library(ggplot2)

#######################
###### Read Data ######
#######################

# Remove problem mesocosms: 44,55,56,63,147,153,166,245,265,473
# Make new cleaned death
# Make nondetected small value biomass
# Merge with div data
# Add block design: 1: endings 1 and 6; 2: 2 and 7; 3: 3 and 8, 4: 4 and 9, 5: 5 and 0

meso_bad <- c(44,55,56,63,147,153,166,245,265,473)

div <- read.csv("data/Ericoid_mesocosm_combinations.csv") %>%
  select(meso_num, p_level, f_level,f_comb) 

dat <- read.table("data/ERM_dataclean_April42025.txt", header = TRUE) %>%
    filter(!meso_num %in% meso_bad) %>%
    mutate(dead = case_when(dead_all == "dead" ~ "dead",
                            dead_all == "dead;no_leaves" ~ "dead",
                            dead_all == "no_plant" ~ "no_plant",
                            dead_all == "sick" ~ "sick",
                            dead_all == "sick;shoot_broke" ~ "sick")) %>%
    mutate(dead = ifelse(is.na(dead), "alive", dead)) %>%
    select(-dead_all) %>%
    #add small value to not detected
    mutate(above_biomass_g = ifelse(above_biomass_g == "not_detected", 0.0001, above_biomass_g)) %>%
    filter(!dead == "dead") %>%
    filter(!dead == "no_plant") %>%
    left_join(div, by = "meso_num") %>%
    mutate(p_level = as.character(p_level), f_level = as.character(f_level), above_biomass_g = as.numeric(above_biomass_g),  belowBEF_biomass_g = as.numeric(belowBEF_biomass_g)) %>%
    mutate(last_digit = meso_num %% 10,
    block = case_when(
      last_digit %in% c(1, 6) ~ 1,
      last_digit %in% c(2, 7) ~ 2,
      last_digit %in% c(3, 8) ~ 3,
      last_digit %in% c(4, 9) ~ 4,
      last_digit %in% c(5, 0) ~ 5)) %>%
  select(-last_digit)

###########################
###### Survival data ######
###########################

dat.all <- read.table("data/ERM_dataclean_April42025.txt", header = TRUE) %>%
  filter(!meso_num %in% meso_bad) %>%
  mutate(dead = case_when(dead_all == "dead" ~ "dead",
                          dead_all == "dead;no_leaves" ~ "dead",
                          dead_all == "no_plant" ~ "no_plant",
                          dead_all == "sick" ~ "sick",
                          dead_all == "sick;shoot_broke" ~ "sick")) %>%
  mutate(dead = ifelse(is.na(dead), "alive", dead)) %>%
  select(-dead_all) %>%
  #add small value to not detected
  mutate(above_biomass_g = ifelse(above_biomass_g == "not_detected", 0.0001, above_biomass_g)) %>%
  #filter(!dead == "dead") %>%
  filter(!dead == "no_plant") %>%
  left_join(div, by = "meso_num") %>%
  mutate(p_level = as.character(p_level), f_level = as.character(f_level), above_biomass_g = as.numeric(above_biomass_g),  belowBEF_biomass_g = as.numeric(belowBEF_biomass_g)) %>%
  mutate(last_digit = meso_num %% 10,
         block = case_when(
           last_digit %in% c(1, 6) ~ 1,
           last_digit %in% c(2, 7) ~ 2,
           last_digit %in% c(3, 8) ~ 3,
           last_digit %in% c(4, 9) ~ 4,
           last_digit %in% c(5, 0) ~ 5)) %>%
  select(-last_digit)

#######################
######## Save #########
#######################

saveRDS(dat, "data/ERM_data_abovegroundind.rds")
saveRDS(dat.all, "data/ERM_data_abovegroundind_all.rds")

# dat a
dat.a <- dat %>%
  select(-belowBEF_biomass_g) %>%
  drop_na(-f_comb)

# ind num per meso
dat.indpermeso <- dat %>% 
  group_by(meso_num) %>%
  summarise(ind_num = n())

#mesocosm level data
meso.a <- dat.a %>%
  group_by(meso_num) %>%
  summarise(aboveBEF_biomass_g = sum(above_biomass_g)) 

dat.meso.a <- meso.a %>%
  left_join(dat.a, by = "meso_num") %>%
  select(meso_num, aboveBEF_biomass_g, p_level, f_level, block) %>%
  distinct(meso_num, .keep_all = TRUE) %>%
  left_join(dat.indpermeso, by = "meso_num")

dat.meso.b <- dat %>%
  select(meso_num, belowBEF_biomass_g, p_level, f_level, block) %>%
  drop_na() %>%
  distinct(meso_num, .keep_all = TRUE) %>%
  left_join(dat.indpermeso, by = "meso_num")

#######################
######## Save #########
#######################

saveRDS(dat.meso.a, "data/ERM_data_abovegroundmeso.rds")
saveRDS(dat.meso.b, "data/ERM_data_belowgroundmeso.rds")

#######################
######## Plot #########
#######################

#### Plot Data ####
#### Ind prod #####

# Summarize means and standard errors
summary_df <- dat.a %>%
  group_by(p_level, f_level) %>%
  summarise(
    mean_biomass = mean(above_biomass_g),
    se = sd(above_biomass_g) / sqrt(n())
  )

plot <- ggplot() +
  geom_jitter(data = dat.a, 
              aes(x = p_level, y = above_biomass_g, color = f_level),
              width = 0.2, alpha = 0.3, size = 1) +
  geom_point(data = summary_df,
             aes(x = p_level, y = mean_biomass, color = f_level),
             position = position_dodge(0.2), size = 2) +
  geom_errorbar(data = summary_df,
                aes(x = p_level, ymin = mean_biomass - se, ymax = mean_biomass + se, color = f_level),
                width = 0.2,
                position = position_dodge(0.2)) +
  scale_color_brewer(palette = "Dark2") +  
  labs(y = "Aboveground Biomass (g)", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  ylim(0, 0.015)  +
  facet_grid(~f_level)

png("figures/ERM_aboveground_ind_data.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

##### Plot Data #####
#### Meso above #####

# Summarize means and standard errors
summary_df <- dat.meso.a %>%
  group_by(p_level, f_level) %>%
  summarise(
    mean_biomass = mean(aboveBEF_biomass_g),
    se = sd(aboveBEF_biomass_g) / sqrt(n())
  )

plot <- 
  ggplot() +
  geom_jitter(data = dat.meso.a, 
              aes(x = p_level, y = aboveBEF_biomass_g, color = f_level),
              width = 0.2, alpha = 0.3, size = 1) +
  geom_point(data = summary_df,
             aes(x = p_level, y = mean_biomass, color = f_level),
             position = position_dodge(0.2), size = 2) +
  geom_errorbar(data = summary_df,
                aes(x = p_level, ymin = mean_biomass - se, ymax = mean_biomass + se, color = f_level),
                width = 0.2,
                position = position_dodge(0.2)) +
  scale_color_brewer(palette = "Dark2") + 
  labs(y = "Aboveground Biomass (g)", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  ylim(0, 0.20)  +
  facet_grid(~f_level)

png("figures/ERM_aboveground_mesocosm_data.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

##### Plot Data #####
#### Meso Below #####

# Summarize means and standard errors
summary_df <- dat.meso.b %>%
  group_by(p_level, f_level) %>%
  summarise(
    mean_biomass = mean(belowBEF_biomass_g),
    se = sd(belowBEF_biomass_g) / sqrt(n())
  )

plot <- 
ggplot() +
  geom_jitter(data = dat.meso.b, 
              aes(x = p_level, y = belowBEF_biomass_g, color = f_level),
              width = 0.2, alpha = 0.3, size = 1) +
  geom_point(data = summary_df,
             aes(x = p_level, y = mean_biomass, color = f_level),
             position = position_dodge(0.2), size = 2) +
  geom_errorbar(data = summary_df,
                aes(x = p_level, ymin = mean_biomass - se, ymax = mean_biomass + se, color = f_level),
                width = 0.2,
                position = position_dodge(0.2)) +
  scale_color_brewer(palette = "Dark2") +  
  labs(y = "Belowground Biomass (g)", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  
  ylim(0, 0.20)  +
  facet_grid(~f_level)

png("figures/ERM_belowground_mesocosm_data.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()
