###########################################################
####### BIOMASS ANALYSES: ABOVE & BELOW, MESO & IND #######
###########################################################

####### PURPOSE OF THIS SCRIPT:
# run biomass models

####### RAW DATA USED:
# ERM_data_abovegroundmeso_predicted.rds
# ERM_data_belowgroundmeso.rds
# ERM_data_abovegroundind_predicted.rds

####### DATA CREATED:
# contrasts_AGB_mesocosm_pdiv.csv
# contrasts_AGB_mesocosm_fdiv.csv
# contrasts_BGB_mesocosm_pdiv.csv
# contrasts_BGB_mesocosm_fdiv.csv
# contrasts_AGB_individual_pdiv.csv
# contrasts_AGB_individual_fdiv.csv
# contrasts_AGB_relyield_pdiv.csv
# contrasts_AGB_relyield_fdiv.csv

###### PLOTS CREATED:
# AGB_mesocosm_pdiv.jpg
# AGB_mesocosm_fdiv.jpg
# BGB_mesocosm_pdiv.jpg
# BGB_mesocosm_fdiv.jpg
# AGB_individual_pdiv.jpg
# AGB_individual_fdiv.jpg
# AGB_pabundance.jpg
# AGB_fabundance.jpg
# AGB_relyield_pdiv.jpg
# AGB_relyield_fdiv.jpg

#######################
#### Load Packages ####
#######################

library(nlstools)
library(tidyverse)
library(lme4)
library(multcompView)
library(lmerTest)
library(emmeans)
library(RColorBrewer)

#######################
###### Read Data ######
#######################

# with predicted data
# only for A as we can correct for each individual based on height
dat.a <- readRDS("data/ERM_data_abovegroundmeso_predicted.rds") %>% rename(AGB = aboveBEF_biomass_g) 
# dat.a <- readRDS("data/ERM_data_abovegroundmeso.rds") %>% rename(AGB = aboveBEF_biomass_g) 

dat.b <- readRDS("data/ERM_data_belowgroundmeso.rds") %>% rename(BGB = belowBEF_biomass_g) 

## set color palettes:
p_colors <- brewer.pal(9, "YlGn")[5:9]
f_colors <- brewer.pal(9, "BuPu")[5:9]

###################
#### MESOCOSM #####
####### AGB #######
###################

# mod
ma <- lmer(log(AGB) ~ p_level*f_level + ind_num + (1|block), dat = dat.a)
anova(ma)

# plot
ref <- emmeans(ma,pairwise~p_level*f_level, data= dat.a, type = "response")
ref.table <- as.data.frame(ref$emmeans) %>% rename(emmean = response)

plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = f_level), method = "glm", se = FALSE) + 
  geom_point(data = dat.a, aes(p_level, AGB, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Aboveground biomass", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) + 
  facet_grid(~f_level) +
  ylim(0,0.10)

png("figures/AGB_mesocosm_pdiv.jpg", width = 10, height = 4, units = 'in', res = 300)
plot
dev.off()

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = p_level), method = "glm", se = FALSE) +
  geom_point(data = dat.a, aes(f_level, AGB, color = p_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Aboveground biomass", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(~p_level) +
  ylim(0,0.10)

png("figures/AGB_mesocosm_fdiv.jpg", width = 10, height = 4, units = 'in', res = 300)
plot
dev.off()

# contrasts
ref <- emmeans(ma, pairwise ~ p_level|f_level, data = dat.a, type = "response")
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_AGB_mesocosm_pdiv.csv", row.names = FALSE)

ref <- emmeans(ma, pairwise ~ f_level|p_level, data = dat.a, type = "response")
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_AGB_mesocosm_fdiv.csv", row.names = FALSE)

###################
#### MESOCOSM #####
####### BGB #######
###################

# mod
mb <- lmer(log(BGB) ~ p_level*f_level + ind_num + (1|block), dat = dat.b)
anova(mb)

# contrasts
ref <- emmeans(mb, pairwise ~ p_level, data = dat.b, type = "response")
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)

# plot
ref <- emmeans(mb,pairwise~p_level*f_level, data= dat.b, type = "response")
ref.table <- as.data.frame(ref$emmeans) %>% rename(emmean = response)

plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = f_level), method = "glm", se = FALSE) +
  geom_point(data = dat.a, aes(p_level, AGB, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Belowground biomass", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) +
  facet_grid(~f_level) +
  ylim(0,0.10)

png("figures/BGB_mesocosm_pdiv.jpg", width = 10, height = 4, units = 'in', res = 300)
plot
dev.off()

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = p_level), method = "glm", se = FALSE) +
  geom_point(data = dat.a, aes(f_level, AGB, color = p_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Belowground biomass", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(~p_level) +
  ylim(0,0.10)

png("figures/BGB_mesocosm_fdiv.jpg", width = 10, height = 4, units = 'in', res = 300)
plot
dev.off()

# contrasts
ref <- emmeans(mb, pairwise ~ p_level, data = dat.b, type = "response")
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_BGB_mesocosm_pdivonly.csv", row.names = FALSE)

ref <- emmeans(mb, pairwise ~ p_level|f_level, data = dat.b, type = "response")
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_BGB_mesocosm_pdiv.csv", row.names = FALSE)

ref <- emmeans(mb, pairwise ~ f_level|p_level, data = dat.b, type = "response")
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_BGB_mesocosm_fdiv.csv", row.names = FALSE)

###################
### IND SPECIES ###
####### AGB #######
###################

dat.a.ind <- readRDS("data/ERM_data_abovegroundind_predicted.rds") %>% select(-AGB) %>% rename(AGB_ind = AGB_pred)

# mod
ma.ind <- lmer(log(AGB_ind) ~ p_level*f_level*plant_sp + (1|block) + (1|block:meso_num), dat = dat.a.ind)
anova(ma.ind)

# plot all 
ref <- emmeans(ma.ind,pairwise~p_level*f_level*plant_sp, data= dat.a.ind, type = "response", pbkrtest.limit = 5000)
ref.table <- as.data.frame(ref$emmeans) %>% rename(emmean = response) 

plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = interaction(plant_sp, f_level)), method = "glm", se = FALSE) +
  #geom_point(data = dat.a.ind, aes(p_level, AGB_ind, color = f_level), 
  #           position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Aboveground biomass", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) +
  facet_grid(plant_sp~f_level, scales = "free_y") 

png("figures/AGB_individual_pdiv.jpg", width = 15, height = 10, units = 'in', res = 300)
plot
dev.off()

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1, position=position_dodge(1)) + 
  geom_smooth(aes(group = interaction(plant_sp, p_level)), method = "glm", se = FALSE) +
  #geom_point(data = dat.a.ind, aes(f_level, AGB_ind, color = p_level), 
  #           position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Aboveground biomass", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(plant_sp~p_level, scales = "free_y") 

png("figures/AGB_individual_fdiv.jpg", width = 15, height = 10, units = 'in', res = 300)
plot
dev.off()

# plot select
ref <- emmeans(ma.ind,pairwise~p_level*f_level*plant_sp, data= dat.a.ind, type = "response", pbkrtest.limit = 5000)
ref.table <- as.data.frame(ref$emmeans) %>% rename(emmean = response) %>%
  filter(plant_sp %in% c("CVU", "GMI"))

plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = interaction(plant_sp, f_level)), method = "glm", se = FALSE) +
  geom_point(data = dat.a.ind %>% filter(plant_sp %in% c("CVU", "GMI")), aes(p_level, AGB_ind, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Aboveground biomass", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 25, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 23),
        axis.text.y=element_text(size = 23),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) +
  facet_grid(plant_sp~f_level, scales = "free_y") +
  ylim(0,0.005)

png("figures/AGB_individual_pdiv_selectspecies.jpg", width = 15, height = 6, units = 'in', res = 300)
plot
dev.off()

ref <- emmeans(ma.ind,pairwise~p_level*f_level*plant_sp, data= dat.a.ind, type = "response", pbkrtest.limit = 5000)
ref.table <- as.data.frame(ref$emmeans) %>% rename(emmean = response) %>%
  filter(plant_sp %in% c("CVU", "VMA"))

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1, position=position_dodge(1)) + 
  geom_smooth(aes(group = interaction(plant_sp, p_level)), method = "glm", se = FALSE) +
  geom_point(data = dat.a.ind %>% filter(plant_sp %in% c("CVU", "VMA")), aes(f_level, AGB_ind, color = p_level), 
            position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Aboveground biomass", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 25, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 23),
        axis.text.y=element_text(size = 23),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(plant_sp~p_level, scales = "free_y") +
  ylim(0,0.02)

png("figures/AGB_individual_fdiv_selectspecies.jpg", width = 15, height = 6, units = 'in', res = 300)
plot
dev.off()

# contrasts
ref <- emmeans(ma.ind, pairwise ~ p_level|plant_sp|f_level, data = dat.a.ind, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_AGB_individual_pdiv.csv", row.names = FALSE)

ref <- emmeans(ma.ind, pairwise ~ f_level|plant_sp|p_level, data = dat.a.ind, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_AGB_individual_fdiv.csv", row.names = FALSE)

###################
### IND PSPECIES ##
###### ABUND ######
###################

# summary
summary <- dat.a.ind %>%
  group_by(p_level, plant_sp, f_level) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(p_level,f_level) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# plot
plot <- ggplot(summary, aes(x = p_level, y = prop, fill = plant_sp)) +
  geom_bar(stat = "identity") +
  labs(x = "Plant diversity", y = "Percent", fill = "Plant species") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "none") +
  scale_fill_viridis_d() + 
  facet_grid(~ f_level)

png("figures/AGB_pabundance.jpg", width = 8, height = 5, units = 'in', res = 300)
plot
dev.off()

###################
### IND FSPECIES ##
###### ABUND ######
###################

# dat
dat.a.ind.f <- dat.a.ind %>%
  filter(!f_level == "0") %>%
  select("plant_sp", "p_level", "f_level", "f_comb") %>%
  separate_rows(f_comb, sep = "_") %>%
  rename(fungal_sp = f_comb)

# summary
summary <- dat.a.ind.f %>%
  group_by(p_level, plant_sp, f_level,fungal_sp) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(p_level,f_level) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# plot
plot <- ggplot(summary, aes(x = p_level, y = prop, fill = fungal_sp)) +
  geom_bar(stat = "identity") +
  labs(x = "Fungal diversity", y = "Percent", fill = "Fungal species") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "none") +
  scale_fill_viridis_d() +
  facet_grid(~ f_level)

png("figures/AGB_fabundance.jpg", width = 8, height = 5, units = 'in', res = 300)
plot
dev.off()

###################
### IND SPECIES ###
## OVERYIELDING ###
###################

# calculate prod in monoculture, p_level = 1
# calculate expected prod in mixture (prop of community)
# compare calculated to read for that species
# first by p div, then p div by f div

# by plant species and f_level
mono.sp.prod.f <- dat.a.ind %>% 
  filter(p_level == "1") %>%
  group_by(meso_num,plant_sp,f_level) %>%
  summarise(prod = sum(AGB_ind)) %>%
  group_by(plant_sp,f_level) %>%
  summarise(mean.prod = mean(prod)) 

# take all !p_level 1 to calculate expected prod
# by plant species and f_level
summary.sp.f <- dat.a.ind %>%
  filter(!p_level == "1") %>%
  group_by(meso_num, p_level, plant_sp, f_level) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(meso_num, p_level,f_level) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  left_join(mono.sp.prod.f, by = c("plant_sp","f_level")) %>%
  mutate(exp.prod = prop*mean.prod) %>%
  mutate(relyield = mean.prod - exp.prod) 

# dat
block <- dat.a %>% select("meso_num", "block")
dat.a.yield <- summary.sp.f %>% left_join(block, by = "meso_num")

# mod
ma.yield <- lmer(relyield ~ p_level*f_level*plant_sp + (1|block) + (1|block:meso_num), dat = dat.a.yield)
anova(ma.yield)

ref <- emmeans(ma.yield,pairwise~p_level*f_level*plant_sp, data= dat.a.yield, type = "response", pbkrtest.limit = 5000)
ref.table <- as.data.frame(ref$emmeans)

# plot
plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =3, alpha = 0.8) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.3,size=1, alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = interaction(plant_sp, f_level)), method = "glm", se = FALSE) +
  geom_point(data = dat.a.yield, aes(p_level, relyield, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Overyielding", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) +
  facet_grid(plant_sp~f_level, scales = "free_y") 

png("figures/AGB_relyield_pdiv.jpg", width = 15, height = 10, units = 'in', res = 300)
plot
dev.off()

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =3, alpha = 0.8) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.3,size=1, alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = interaction(plant_sp, p_level)), method = "glm", se = FALSE) +
  geom_point(data = dat.a.yield, aes(f_level, relyield, color = p_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Overyielding", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(plant_sp~p_level, scales = "free_y") 

png("figures/AGB_relyield_fdiv.jpg", width = 15, height = 10, units = 'in', res = 300)
plot
dev.off()

# contrasts
ref <- emmeans(ma.yield, pairwise ~ p_level|plant_sp|f_level, data = dat.a.yield, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_AGB_relyield_pdiv.csv", row.names = FALSE)

ref <- emmeans(ma.yield, pairwise ~ f_level|plant_sp|p_level, data = dat.a.yield, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_AGB_relyield_fdiv.csv", row.names = FALSE)
