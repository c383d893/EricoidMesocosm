######################################################
###### STATS AND PLOTS ON MGR INCLUDING F_LEVEL ###### 
######################################################

####### PURPOSE OF THIS SCRIPT:
# model MGR using an LMM including fungal_species*plant_sp*p_level
# extract emmeans and contrasts to compare MGR between variables of interest

# model MGR using an LMM including f_level*plant_sp*p_level
# extract emmeans and contrasts to compare MGR between variables of interest

###### DATA USED:
# MGR_biomass_prep.rds
# MGR_biomass_incl_fspecies_prep.rds

###### DATA CREATED:
# contrasts_MGR_plant_sp_x_p_level.csv --> contains contrasts of MGR between p_level of one plant_sp
# contrasts_MGR_3way_flevel.csv --> contains contrasts of MGR between f_levels for each p_level of one plant_sp
# contrasts_MGR_3way_plevel.csv --> contains contrasts of MGR between p_levels for each f_level of one plant_sp

###### PLOTS CREATED: 
# MGR_per_plevel_per_sp.jpg --> plot displaying MGR per p_level per plant_sp; individual y-axis for each small plot
# Masterplot_MGR_flevel.jpg --> plot displaying MGR for each individual combination of f_level*plant_sp*p_level; f_level on x-axis
# Masterplot_MGR_plevel.jpg --> plot displaying MGR for each individual combination of f_level*plant_sp*p_level; p_level on x-axis

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

###### load data 
biomass_df <- readRDS("data/MGR_biomass_incl_fspecies_prep.rds")

######  only select non-sterile samples 
biomass_df_nstr = subset(biomass_df, biomass_df[,"f_level"]!= 0)

######  convert columns into factors (so it's not a continuous variable where the mean is calculated afterwards)
biomass_df_nstr$p_level <- as.factor(biomass_df_nstr$p_level)
biomass_df_nstr$block <- as.factor(biomass_df_nstr$block)
biomass_df_nstr$meso_num <- as.factor(biomass_df_nstr$meso_num)
biomass_df_nstr$plant_sp <- as.factor(biomass_df_nstr$plant_sp)
biomass_df_nstr$f_level <- as.factor(biomass_df_nstr$f_level)
biomass_df_nstr$plant_id <- as.factor(biomass_df_nstr$plant_id)
biomass_df_nstr$f_comb <- as.factor(biomass_df_nstr$plant_id)

#################################
###### feed data into model #####
#################################

######  get smallest value in dataset
min_value <- min(biomass_df_nstr$response)
min_value

biomass_df_nstr <- biomass_df_nstr %>% mutate(response = response + 1)

###### model 1: plant_sp*plant_div*fung_div
set.seed(123)
m_psp_plevel_fsp = lmer(log(biomass_df_nstr$response) ~ f_comb*plant_sp*p_level + height_i_cm + (1|block) + (1|block:meso_num), data=biomass_df_nstr)      
summary(m_psp_plevel_fsp)
anova(m_psp_plevel_fsp)

###### model 2: plant_sp*p_level + f_comb*p_level + f_comb*plant_sp
set.seed(123)
m_psp_plevel_fsp_2 = lmer(log(biomass_df_nstr$response) ~ plant_sp*p_level + f_comb*p_level + f_comb*plant_sp + height_i_cm + (1|block) + (1|block:meso_num), data=biomass_df_nstr)      
summary(m_psp_plevel_fsp_2)
anova(m_psp_plevel_fsp_2)

###### model 3: plant_sp*p_level + f_comb
set.seed(123)
m_psp_plevel_fsp_3 = lmer(log(biomass_df_nstr$response) ~ f_comb + plant_sp*p_level + height_i_cm + (1|block) + (1|block:meso_num), data=biomass_df_nstr)      
summary(m_psp_plevel_fsp_3)
anova(m_psp_plevel_fsp_3)

#########################################
###### extract emmeans & contrasts ######
#########################################

##########################
###### plant_sp*plant_div 

###### extract emmeans
means_total = emmeans::emmeans(m_psp_plevel_fsp_3, ~ plant_sp*p_level, type = "response")   

###### get summary of emmeans
summary_means_total <- summary(means_total)
summary_means_total <- na.omit(summary_means_total) 

###### extract contrasts between plant diversity levels per plant species
ref <- emmeans(m_psp_plevel_fsp_3, pairwise ~ p_level|plant_sp, data = biomass_df_nstr, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df , file = "results/contrasts_MGR_plant_sp_x_p_level.csv", row.names = FALSE)

#############################################
##### plot MGR per p_level per plant_sp #####
#############################################

###### get compact letter display
letters_df <- cld(means_total, 
                  by = "plant_sp", 
                  adjust = "sidak", 
                  Letters = letters) %>%
  as.data.frame() %>%
  dplyr::select(plant_sp, p_level, .group)


###### convert means to data frame for plotting
emm_df <- as.data.frame(means_total) %>%
  transmute(
    plant_sp,
    p_level,
    response,      
    SE         = SE,
    df         = df,
    lower.CL   = lower.CL,
    upper.CL   = upper.CL
  ) %>%
  left_join(letters_df, by = c("plant_sp", "p_level"))


#########################################
###### plot with individual y-axis ranges

emm_df$p_level <- factor(emm_df$p_level, levels = sort(unique(emm_df$p_level)))

emm_df <- emm_df %>%
  group_by(plant_sp) %>%
  mutate(label_y = upper.CL + 0.05 * (max(response, na.rm = TRUE) - min(response, na.rm = TRUE))) %>%
  ungroup()

plot <- ggplot(emm_df, aes(x = p_level, y = response)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = .group, y = label_y), vjust = 0) +
  facet_wrap(~ plant_sp, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  
  labs(
    x = "Plant diversity",
    y = "MGR",
    title = "EMMs of MGR by plant diversity level within each species",
    size = 12
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

png("figures/MGR_per_plevel_per_sp.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

############################
###### Fungal level ########
############################

#######################
###### Read Data ######
#######################

###### clear workspace
rm(list = ls())

###### load data 
biomass_df <- readRDS("data/MGR_biomass_prep.rds")

######  only select non-sterile samples 
biomass_df_nstr = subset(biomass_df, biomass_df[,"f_level"]!= 0)

######  convert columns into factors 
biomass_df_nstr$p_level <- as.factor(biomass_df_nstr$p_level)
biomass_df_nstr$block <- as.factor(biomass_df_nstr$block)
biomass_df_nstr$meso_num <- as.factor(biomass_df_nstr$meso_num)
biomass_df_nstr$plant_sp <- as.factor(biomass_df_nstr$plant_sp)
biomass_df_nstr$f_level <- as.factor(biomass_df_nstr$f_level)
biomass_df_nstr$plant_id <- as.factor(biomass_df_nstr$plant_id)

#################################
###### feed data into model #####
#################################

######  get smalles value in dataset
min_value <- min(biomass_df_nstr$response, na.rm = TRUE)
min_value

biomass_df_nstr <- biomass_df_nstr %>% mutate(response = response + 1)

###### model: plant_sp*plant_div*fung_div
set.seed(123)
m_psp_plevel_flevel = lmer(log(biomass_df_nstr$response) ~ f_level*plant_sp*p_level + height_i_cm + (1|block) + (1|block:meso_num), data=biomass_df_nstr)      
summary(m_psp_plevel_flevel)
anova(m_psp_plevel_flevel)

#########################################
###### extract emmeans & contrasts ######
#########################################

################################
#### fung_div*plant_sp*plant_div

# extract emmeans 
means_total <- emmeans::emmeans(
  m_psp_plevel_flevel, 
  ~ f_level * plant_sp * p_level, 
  type = "response", 
  pbkrtest.limit = 4000, 
  lmerTest.limit = 4000
)

# make summary of emmeans
summary_means_total <- summary(means_total)
summary_means_total <- na.omit(summary_means_total) # excludes PJA p_level 2

#################################
#### f_level | plant_sp | p_level

# extract emmeans
means_sp_per_flevel <- emmeans::emmeans(
  m_psp_plevel_flevel,
  ~ f_level | plant_sp * p_level,  
  type = "response",
  pbkrtest.limit = 4000, 
  lmerTest.limit = 4000
)

# make summary of emmeans
summary_means_sp_per_flevel <- summary(means_sp_per_flevel)
summary_means_sp_per_flevel <- na.omit(summary_means_sp_per_flevel) # excludes PJA p_level 2

# extract contrasts between f_level within each p_level of one plant_sp
ref <- emmeans(m_psp_plevel_flevel, pairwise ~ f_level | plant_sp * p_level, data = biomass_df_nstr, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_MGR_3way_flevel.csv", row.names = FALSE)

#################################
#### p_level | plant_sp | f_level

# extract emmeans 
means_sp_per_plevel <- emmeans::emmeans(
  m_psp_plevel_flevel,
  ~ p_level | plant_sp * f_level,  
  type = "response",
  pbkrtest.limit = 4000, 
  lmerTest.limit = 4000
)

# make summary of emmeans
summary_means_sp_per_plevel <- summary(means_sp_per_plevel)
summary_means_sp_per_plevel <- na.omit(summary_means_sp_per_plevel) # excludes PJA p_level 2

# extract contrasts between p_level within each f_level of one plant_sp
ref <- emmeans(m_psp_plevel_flevel, pairwise ~ p_level | plant_sp * f_level, data = biomass_df_nstr, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_MGR_3way_plevel.csv", row.names = FALSE)

#######################################################
########## PLOT BMGR MASTERPLOT PER F_LEVEL ###########
#######################################################

# define a consistent color palette for fung_sp
f_level_colors <- c(
  "1" = "#fc9403",
  "2" = "#379e20",
  "4" = "#7f209e",
  "7" = "#03d7fc"
)


# add a column indicating if CI crosses zero
summary_means_sp_per_flevel <- summary_means_sp_per_flevel %>%
  mutate(cl_crosses_zero = ifelse(lower.CL < 0 & upper.CL > 0, "dashed", "solid"))

# create a variable for line types
summary_means_sp_per_flevel$line_type <- ifelse(
  summary_means_sp_per_flevel$lower.CL < 0 & summary_means_sp_per_flevel$upper.CL > 0, 
  "yes", 
  "no"
)

# create letters per facet
cld_facet <- cld(
  means_sp_per_flevel,               
  by       = c("plant_sp", "p_level"),      # do a separate CLD in every facet
  Letters  = letters,                       
  adjust   = "sidak"
) %>% 
  as_tibble() %>%                           
  mutate(
    .group = str_trim(.group),              
    f_level = factor(f_level, levels = c("1","2","4","7"))
  )

# merge with the summary that feeds ggplot 
summary_means_sp_per_flevel_letters <- summary_means_sp_per_flevel %>% 
  left_join(
    cld_facet %>% 
      dplyr::select(f_level, plant_sp, p_level, .group),
    by = c("f_level", "plant_sp", "p_level")
  ) %>%
  group_by(plant_sp, p_level) %>%
  mutate(
    # count number of unique letters in this facet
    unique_groups = n_distinct(.group)
  ) %>%
  ungroup() %>%
  # keep only rows where >1 letter group exists in that facet
  filter(unique_groups > 1)

# make plot 
plot <- ggplot(summary_means_sp_per_flevel, aes(x = f_level, y = response, fill = f_level, color = f_level)) +
  geom_text(
    data = summary_means_sp_per_flevel_letters,
    aes(x = f_level, label = .group),
    y = 8.5,               # fixed height above your y-limit of 5
    color = "black",
    size       = 3,
    vjust      = 0,          # text sits on top of y_lab
    show.legend = FALSE
  )+
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  geom_pointrange(
    aes(
      ymin = lower.CL, ymax = upper.CL, 
      linetype = line_type  # Map the line type
    ),
    fatten = 3.5, 
    size = 0.6, 
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_manual(values = f_level_colors) +
  scale_color_manual(values = f_level_colors) +
  scale_linetype_manual(
    values = c("yes" = "dashed", "no" = "solid"),  # Map line types
    guide = guide_legend(override.aes = list(
      color = "black", 
      size = 0.5, 
      fill = NA, 
      shape = NA  # Remove points from the legend
    ))
  ) +
  scale_y_continuous(
    limits = c(-3, 9),
    sec.axis = sec_axis(transform = ~., name = "Plant diversity level", breaks = NULL)
  ) +
  facet_grid(p_level ~ plant_sp, scales = "fixed") +
  geom_point(
    data = biomass_df_nstr, 
    aes(x = f_level, y = response, fill = f_level, color = f_level), 
    alpha = 0.2, 
    size = 0.9,
    position = position_jitter(width = 0.2, height = 0) 
  ) +
  labs(
    title = "MGR",
    x = "Fungal diversity level",
    y = "MGR",
    linetype = "CI crosses 0"  # Legend title for line type
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(, hjust = 1, size = 12),
    strip.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 30)
  )

png("figures/Masterplot_MGR_flevel.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()

#######################################################
########## PLOT BMGR MASTERPLOT PER P_LEVEL ###########
#######################################################

# define a consistent color palette for fung_sp
p_level_colors <- c(
  "1" = "#fc9403",
  "2" = "#379e20",
  "4" = "#7f209e",
  "8" = "#03d7fc"
)


# add a column indicating if CI crosses zero
summary_means_sp_per_plevel <- summary_means_sp_per_plevel %>%
  mutate(cl_crosses_zero = ifelse(lower.CL < 0 & upper.CL > 0, "dashed", "solid"))

# create a variable for line types
summary_means_sp_per_plevel$line_type <- ifelse(
  summary_means_sp_per_flevel$lower.CL < 0 & summary_means_sp_per_plevel$upper.CL > 0, 
  "yes", 
  "no"
)

# create letters per facet
cld_facet_p <- cld(
  means_sp_per_plevel,               
  by       = c("plant_sp", "f_level"),      # do a separate CLD in every facet
  Letters  = letters,                       
  adjust   = "sidak"
) %>% 
  as_tibble() %>%                           
  mutate(
    .group = str_trim(.group),             
    p_level = factor(p_level, levels = c("1","2","4","8"))
  )

# merge with the summary that feeds ggplot 
summary_means_sp_per_plevel_letters <- summary_means_sp_per_plevel %>% 
  left_join(
    cld_facet_p %>% 
      dplyr::select(p_level, plant_sp, f_level, .group),
    by = c("p_level", "plant_sp", "f_level")
  ) %>%
  group_by(plant_sp, f_level) %>%
  mutate(
    # count number of unique letters in this facet
    unique_groups = n_distinct(.group)
  ) %>%
  ungroup() %>%
  # keep only rows where >1 letter group exists in that facet
  filter(unique_groups > 1)

# make plot 
plot <- ggplot(summary_means_sp_per_plevel, aes(x = p_level, y = response, fill = p_level, color = p_level)) +
  geom_text(
    data = summary_means_sp_per_plevel_letters,
    aes(x = p_level, label = .group),
    y = 8.5,               # fixed height above your y-limit of 5
    color = "black",
    size       = 3,
    vjust      = 0,          # text sits on top of y_lab
    show.legend = FALSE
  )+
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  geom_pointrange(
    aes(
      ymin = lower.CL, ymax = upper.CL, 
      linetype = line_type  # Map the line type
    ),
    fatten = 3.5, 
    size = 0.6, 
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_manual(values = p_level_colors) +
  scale_color_manual(values = p_level_colors) +
  scale_linetype_manual(
    values = c("yes" = "dashed", "no" = "solid"),  # Map line types
    guide = guide_legend(override.aes = list(
      color = "black", 
      size = 0.5, 
      fill = NA, 
      shape = NA  # Remove points from the legend
    ))
  ) +
  scale_y_continuous(
    limits = c(-3, 9),
    sec.axis = sec_axis(transform = ~., name = "Fungal diversity level", breaks = NULL)
  ) +
  facet_grid(f_level ~ plant_sp, scales = "fixed") +
  geom_point(
    data = biomass_df_nstr, 
    aes(x = p_level, y = response, fill = p_level, color = p_level), 
    alpha = 0.2, 
    size = 0.9,
    position = position_jitter(width = 0.2, height = 0) 
  ) +
  labs(
    title = "MGR",
    x = "Plant diversity level",
    y = "MGR",
    linetype = "CI crosses 0"  # Legend title for line type
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(, hjust = 1, size = 12),
    strip.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 30)
  )

png("figures/Masterplot_MGR_plevel.jpg", width = 10, height = 10, units = 'in', res = 300)
plot
dev.off()
