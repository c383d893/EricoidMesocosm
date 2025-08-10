######################################################
#### STATS AND PLOTS COMPLEMENTARITY V. SELECTION #### 
######################################################

####### PURPOSE OF THIS SCRIPT:
# calculate complementarity and selection
# determine whether this varies based on p_level or f_level

###### DATA USED:
# ERM_data_abovegroundind_predicted.rds

###### DATA CREATED:
# PLOT.BM.g.m2.2025SECENE.csv --> contains complementarity and selection values per species and f_level
# contrasts_CE_AGB_individual_pdiv.csv --> 
# contrasts_CE_AGB_individual_fdiv.csv --> 
# contrasts_SE_AGB_individual_pdiv.csv --> 
# contrasts_SE_AGB_individual_fdiv.csv --> 

###### PLOTS CREATED: 
# CE_AGB_individual_pdiv.jpg --> CE pdiv
# CE_AGB_individual_fdiv.jpg --> CE fdiv
# SE_AGB_individual_pdiv.jpg --> SE pdiv
# SE_AGB_individual_fdiv.jpg --> SE fdiv

#######################
#### Load Packages ####
#######################

library(tidyverse)
library(emmeans)
library(RColorBrewer)
source("functions/SECENEfnct.R")

# set color palettes:
p_colors <- brewer.pal(9, "YlGn")[5:9]
f_colors <- brewer.pal(9, "BuPu")[5:9]

#######################
###### Read Data ######
#######################

Biomass2025 <- readRDS("data/ERM_data_abovegroundind_predicted.rds")
Biomass2025$plant_spf_level <- paste(Biomass2025$plant_sp,Biomass2025$f_level,sep="")
Biomass2025$plant_spf_level <- as.factor(Biomass2025$plant_spf_level)

#2023 Calculating Mi, Yoj, Yo, and RYei

#M[i] is the yield of speces i in monoculture 
# subset to monoculture
M.idf <- subset(Biomass2025,p_level==1)
# aggregate biomass by species and treatment, rename cols, make sp by treat col
M.idfmean <- aggregate(M.idf$AGB,by=list(M.idf$plant_sp,M.idf$f_level),mean,na.rm=T)
names(M.idfmean)[1:3] <- paste(c("plant_sp","f_level","M.i.g"))
M.idfmean$plant_spf_level <- paste(M.idfmean$plant_sp,M.idfmean$f_level,sep="")
M.idfmean$plant_spf_level <- as.factor(M.idfmean$plant_spf_level)
M.idfmean <- M.idfmean[,-c(1:2)]
# join with biomass data by sp by treat col
Biomass2025 <- left_join(Biomass2025,M.idfmean,by="plant_spf_level")

#Yoj is the total observed yield each species in the mixture, which is coded as SppBiomass.g.m2

#Yo is the total observed yield of the mixture
# aggregate biomass by plot no, rename cols
TOTMESOCOVER <- aggregate(Biomass2025$AGB,by=list(Biomass2025$meso_num),sum)
names(TOTMESOCOVER)[1:2]<-paste(c("meso_num","BMTOTMESO.g"))
# join with biomass data by sp by treat col
Biomass2025 <- left_join(Biomass2025,TOTMESOCOVER, by = "meso_num")

#RYej is the expected relative yield of the species in mixture, which is the proportion seeded/planted
Biomass2025$p_level <- as.numeric(Biomass2025$p_level)
Biomass2025$PROPCOMP<-(1/Biomass2025$p_level)
Biomass2025$PROPCOMP100<-Biomass2025$PROPCOMP*100

#Running the CE, SE, and NE on the dataframe
Biomass2025$SPPNO <- Biomass2025$p_level
BM.g.2025SECENE<-SECE(df=Biomass2025,
                      PLOTNO="meso_num",
                      Mi="M.i.g",
                      Yoi="AGB",
                      Yo="BMTOTMESO.g",
                      RYej="PROPCOMP")

#This step is very important. Biodiversity effects describe mixed communities, not monocultures. 
#I have found it easier to do the SE, CE, and NE calculations with the monocultures in the dataset, so I can more easily calculate Mi.
#However, the tradeoff is that we have to remember to remove the monocultures from the dataset ASAP, because the calculations yield values but they are not valid.
#So, it's a good idea to remove the monocultures from the dataset before we save our calculations
BM.g.2025SECENE.mix <- subset(BM.g.2025SECENE,SPPNO!=1&SPPNO!=0)

#removing duplicated values for total plot biomass, NE (deltaY), CE, and SE for analysis
PLOT.BM.g.2025SECENE <- BM.g.2025SECENE.mix[,c(match(c("meso_num","BMTOTMESO.g","deltaY","CE","SE"),names(BM.g.2025SECENE.mix)))]
PLOT.BM.g.2025SECENE <- PLOT.BM.g.2025SECENE[!duplicated(PLOT.BM.g.2025SECENE), ]

#save
write.csv(PLOT.BM.g.2025SECENE, "data/PLOT.BM.g.m2.2025SECENE.csv")

cs.dat <- Biomass2025 %>% select(meso_num, p_level, f_level,plant_sp, block) %>% left_join(PLOT.BM.g.2025SECENE,by="meso_num") %>% filter(!p_level == 1) %>% mutate(p_level = as.factor(p_level))

###################
##### CE MODEL ####
###################

mc <- lmer(CE ~ p_level*f_level + (1|block), dat = cs.dat)
anova(mc)

# plot
ref <- emmeans(mc,pairwise~p_level*f_level, data= cs.dat, type = "response", pbkrtest.limit = 5000)
ref.table <- as.data.frame(ref$emmeans) 

plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = f_level), method = "glm", se = FALSE) + 
  geom_point(data = cs.dat, aes(p_level, CE, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Complementarity", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) + 
  facet_grid(~f_level) +
  ylim(0,0.30) 

png("figures/CE_AGB_individual_pdiv.jpg", width = 15, height = 6, units = 'in', res = 300)
plot
dev.off()

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = p_level), method = "glm", se = FALSE) +
  geom_point(data = cs.dat, aes(f_level, CE, color = p_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Complementarity", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(~p_level) +
  ylim(0,0.10)

png("figures/CE_AGB_individual_fdiv.jpg", width = 15, height = 6, units = 'in', res = 300)
plot
dev.off()

# contrasts
ref <- emmeans(mc, pairwise ~ p_level|f_level, data = cs.dat, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_CE_AGB_individual_pdiv.csv", row.names = FALSE)

ref <- emmeans(mc, pairwise ~ f_level|p_level, data = cs.dat, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_CE_AGB_individual_fdiv.csv", row.names = FALSE)

###################
##### SE MODEL ####
###################

ms <- lmer(SE ~ p_level*f_level + (1|block), dat = cs.dat)
anova(ms)

# plot
ref <- emmeans(ms,pairwise~p_level*f_level, data= cs.dat, type = "response", pbkrtest.limit = 5000)
ref.table <- as.data.frame(ref$emmeans) 

plot <- ggplot(ref.table, aes(p_level, emmean, color=f_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = f_level), method = "glm", se = FALSE) + 
  geom_point(data = cs.dat, aes(p_level, SE, color = f_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Selection", x = "Plant diversity", color = "Fungal diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = f_colors) + 
  facet_grid(~f_level) +
  ylim(-0.5,0.5) 

png("figures/SE_AGB_individual_pdiv.jpg", width = 15, height = 6, units = 'in', res = 300)
plot
dev.off()

plot <- ggplot(ref.table, aes(f_level, emmean, color=p_level)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.3,size=1,alpha = 0.8, position=position_dodge(1)) + 
  geom_smooth(aes(group = p_level), method = "glm", se = FALSE) +
  geom_point(data = cs.dat, aes(f_level, SE, color = p_level), 
             position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw()+
  scale_y_continuous()+
  labs(y = "Selection", x = "Fungal diversity", color = "Plant diversity") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        axis.text.y=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = p_colors) +
  facet_grid(~p_level) +
  ylim(-0.5,0.5)

png("figures/SE_AGB_individual_fdiv.jpg", width = 15, height = 6, units = 'in', res = 300)
plot
dev.off()

# contrasts
ref <- emmeans(ms, pairwise ~ p_level|f_level, data = cs.dat, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_SE_AGB_individual_pdiv.csv", row.names = FALSE)

ref <- emmeans(ms, pairwise ~ f_level|p_level, data = cs.dat, type = "response", pbkrtest.limit = 5000)
pairwise_df <- as.data.frame(ref$contrasts) %>% filter(p.value < 0.05)
write.csv(pairwise_df, file = "results/contrasts_SE_AGB_individual_fdiv.csv", row.names = FALSE)
