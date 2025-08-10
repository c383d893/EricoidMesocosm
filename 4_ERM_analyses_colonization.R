###########################################################
####### MYC COLONIZATION ASSESSMENT: BY FUNGAL LEVEL ######
###########################################################

####### PURPOSE OF THIS SCRIPT:
# test myc colonization between f-levels

####### RAW DATA USED:
# ERM_colonization2025.txt

####### DATA CREATED:
# NA

###### PLOTS CREATED:
# colonization_sterilevs.inoc.jpg
# colonization_sterilevs.f_level

#######################
#### Load Packages ####
#######################

library(tidyverse)
library(emmeans)
library(RColorBrewer)

## set color palettes:
f_colors <- brewer.pal(9, "BuPu")[5:9]

#######################
###### Read Data ######
#######################

coldat <- read.table("data/ERM_colonization2025.txt", header = TRUE) %>% 
  drop_na() %>%
  mutate(f_level = as.character(f_level)) %>%
  mutate(f_level2 = as.character(ifelse(f_level == 0, 0, 1)))

#######################
#### Model/Plot #######
#######################

# sterile v inoculated
colmod <- glm(cbind(colonized,not) ~ f_level2, data = coldat, family = "binomial")
summary(colmod)
anova(colmod)

# plot
ref <- emmeans(colmod, pairwise~f_level2, data = coldat, type = "response")
ref.table <- as.data.frame(ref$emmeans) 

plot <- ggplot(ref.table, aes(f_level2, prob, color=f_level2)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_point(data = coldat, aes(f_level2, percent, color = f_level2), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  scale_y_continuous()+
  labs(y = "Root colonization (%)", x = "Fungal diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "none") +
  scale_color_manual(values = f_colors) +
  ylim(0,0.8)

png("figures/colonization_sterilevs.inoc.jpg", width = 6, height = 5, units = 'in', res = 300)
plot
dev.off()

# sterile v inoculated by level
colmod <- glm(cbind(colonized,not) ~ f_level, data = coldat, family = "binomial")
summary(colmod)
anova(colmod)

# plot
ref <- emmeans(colmod, pairwise~f_level, data = coldat, type = "response")
ref.table <- as.data.frame(ref$emmeans)

plot <- ggplot(ref.table, aes(f_level, prob, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_point(data = coldat, aes(f_level, percent, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Root colonization (%)", x = "Fungal diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "none") +
  scale_color_manual(values = f_colors) +
  ylim(0,1.0)

png("figures/colonization_sterilevs.f_level.jpg", width = 6, height = 5, units = 'in', res = 300)
plot
dev.off()
