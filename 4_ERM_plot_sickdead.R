###################################################
###### STATS AND PLOTS ON SICK / DEAD PLANTS ###### 
###################################################

###### PURPOSE OF THIS SCRIPT:
# generate several plots to get an overview of how many plants per plant species and per plant or fungal diversity were sick / dead

###### DATA USED:
# ERM_data_abovegroundind_all.rds

###### PLOTS CREATED: 
# Heatmap_dead_flevel.png --> heatmap displaying the number of dead plants for each combination of plant species x fungal diversity level; f_level on x-axis
# Heatmap_dead_fspecies.png --> heatmap displaying the number of dead plants for each combination of plant species x fungal species; f_level on x-axis; fungal species on x-axis

##############################
###### install packages ######
##############################

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggpattern)

#################################
###### load & prepare data ######
#################################

###### load data 
ERM_data <- readRDS("data/ERM_data_abovegroundind_all.rds")

######  write plant & fungal species in full length 
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="PJA", replacement="Pieris japonica", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="GPR", replacement="Gaultheria procumbens", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="VCO", replacement="Vaccinium corymbosum", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="VVI", replacement="Vaccinium vitis-idaea", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="KLA", replacement="Kalmia latifolia", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="GSH", replacement="Gaultheria shallon", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="VAN", replacement="Vaccinium angustifolium", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="VMY", replacement="Vaccinium myrtillus", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="CVU", replacement="Calluna vulgaris", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="VMA", replacement="Vaccinium macrocarpon", plant_sp))
ERM_data = transform(ERM_data, plant_sp=gsub(pattern="GMI", replacement="Gaultheria miqueliana", plant_sp))

ERM_data = transform(ERM_data, f_comb=gsub(pattern="JPK132", replacement="Serendipitaceae sp.", f_comb))
ERM_data = transform(ERM_data, f_comb=gsub(pattern="MEB", replacement="Hyaloscypha bicolor", f_comb))
ERM_data = transform(ERM_data, f_comb=gsub(pattern="MEV", replacement="Hyaloscypha variabilis", f_comb))
ERM_data = transform(ERM_data, f_comb=gsub(pattern="MGR", replacement="Hyaloscypha gryndleri", f_comb))
ERM_data = transform(ERM_data, f_comb=gsub(pattern="OMA", replacement="Oidiodendron maius", f_comb))
ERM_data = transform(ERM_data, f_comb=gsub(pattern="RER", replacement="Hyaloscypha hepaticicola_1", f_comb))
ERM_data = transform(ERM_data, f_comb=gsub(pattern="RHE", replacement="Hyaloscypha hepaticicola_2", f_comb))

####### remove unnecessary columns 
ERM_data <- ERM_data %>% select(-c("height_i_cm", "height_f_cm", "above_biomass_g", "belowBEF_biomass_g"))

###### count alive plants
total_alive <- ERM_data[(ERM_data$dead %in% c("alive")), ]
# 4596

###### count sick plants
total_sick <- ERM_data[(ERM_data$dead %in% c("sick")), ]
# 94

###### count dead plants
total_dead <- ERM_data[(ERM_data$dead %in% c("dead")), ]
# 111

###### define order of plant species within plots
plant_order <- c(
  "Calluna vulgaris",
  "Gaultheria miqueliana",
  "Gaultheria procumbens",
  "Gaultheria shallon",
  "Kalmia latifolia",
  "Pieris japonica",
  "Vaccinium angustifolium",
  "Vaccinium corymbosum",
  "Vaccinium macrocarpon",
  "Vaccinium myrtillus",
  "Vaccinium vitis-idaea"
)


##########################################################################
###### heatmap of dead plants per species per fungal diversity level #####
##########################################################################

###### count number of dead plants per species per f_level
df_dead_flevel <- total_dead %>%
  count(plant_sp, f_level, name = "n_dead")

###### add missing plant species (they are not in the total_dead because none of the plants of these species died)
# ensure all unique p_levels are included
all_flevels <- df_dead_flevel %>% distinct(f_level) %>% pull(f_level)

# create a new data frame for Pieris japonica across all p_levels
add_rows_2 <- expand.grid(
  plant_sp = c("Pieris japonica", "Vaccinium macrocarpon", "Vaccinium myrtillus"),
  f_level = all_flevels,
  n_dead = NA_real_
) %>%
  as_tibble()

# add to main data frame
df_dead_flevel_expanded <- bind_rows(df_dead_flevel, add_rows_2)

###### fill in 0 for missing values
df_heatmap_dead_flevel <- df_dead_flevel_expanded %>%
  complete(plant_sp, f_level, fill = list(n_dead = 0)) %>%
  mutate(
    is_zero = n_dead == 0,
    is_na = FALSE          
  )

###### sort heatmap according to specified order of plant species
df_heatmap_dead_flevel$plant_sp <- factor(df_heatmap_dead_flevel$plant_sp, levels = plant_order)

###### plot
heatmap_dead_flevel <- ggplot(df_heatmap_dead_flevel, aes(x = f_level, y = plant_sp)) +
  geom_tile(aes(fill = n_dead), color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "black",
    name = "Dead Plants"
  ) +
  labs(
    x = "Fungal Diversity",
    y = "Plant Species",
    title = "Heatmap of Dead Plants by Fungal Diversity"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

png("figures/Heatmap_dead_flevel.jpg", width = 6.5, height = 6, units = 'in', res = 300)
heatmap_dead_flevel
dev.off()


########################################################################
###### heatmap of dead plants per plant species per fungal species #####
########################################################################

###### extract all rows with fungal diversity level 1
ERM_data_flevel_1 <- ERM_data %>%
  filter(f_level == 1)

unique_combos_1 <- ERM_data_flevel_1 %>% distinct(plant_sp, f_comb)
# missing combination: PJA-OMA

###### count number of dead plants per plant species per fungal species within f_level 1
df_dead_fspecies_flevel1 <- total_dead %>%
  filter(f_level == 1) %>%
  count(plant_sp, f_comb, name = "n_dead")

###### fill in all missing plant_sp-f_species combinations and add n_dead = 0
df_dead_fspecies_flevel1_completed <- df_dead_fspecies_flevel1 %>%
  complete(
    plant_sp = unique(ERM_data_flevel_1$plant_sp),
    f_comb   = unique(ERM_data_flevel_1$f_comb),
    fill = list(n_dead = 0)
  )

df_dead_fspecies_flevel1_completed$n_dead[
  df_dead_fspecies_flevel1_completed$plant_sp == "Pieris japonica" & 
    df_dead_fspecies_flevel1_completed$f_comb == "Oidiodendron maius"
] <- NA_real_

###### sort heatmap according to specified order of plant species
df_dead_fspecies_flevel1_completed$plant_sp <- factor(df_dead_fspecies_flevel1_completed$plant_sp, levels = plant_order)

###### plot
heatmap_dead_fspecies <- ggplot(df_dead_fspecies_flevel1_completed, aes(x = f_comb, y = plant_sp)) +
  # Main heat map for non-NA tiles
  geom_tile(aes(fill = n_dead), color = "black", data = subset(df_dead_fspecies_flevel1_completed, !is.na(n_dead))) +
  # Hatched overlay for n_dead == NA; e.g. no data present for this specific plant-fungus combination
  geom_tile_pattern(
    data = subset(df_dead_fspecies_flevel1_completed, is.na(n_dead)),
    pattern = "stripe",
    fill = "white",
    color = "black",
    pattern_size = 0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    name = "Dead Plants",
  ) +
  labs(
    x = "Fungal Species",
    y = "Plant Species",
    title = "Heatmap of Dead Plants by Fungal Species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("figures/Heatmap_dead_fspecies.jpg", width = 6.5, height = 6, units = 'in', res = 300)
heatmap_dead_fspecies
dev.off()
