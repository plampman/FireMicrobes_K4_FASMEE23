#########################################################################
###
### Program name: FASMEE23 and Konza4 comparison/combined models
###
### Purpose: FASMEE23 and Konza4 mixed models
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/01/2025
###
#########################################################################

library(tidyverse)
library(emmeans)
library(lme4)
library(nlme)
library(glmmTMB)
library(DHARMa)
library(mgcv)
# library(vegan)
library(ggstatsplot)

#Spores
#------------------------------------------------------------------------------------------

k4_spores <- read.csv('Konza4_spores_PA_C_20250501.csv', header = T)
fasmee23_spores <- read.csv('FASMEE23_Spores_PA_C_20250501.csv', header = T) #Blue only

k4_spores <- k4_spores %>%
  select(Project, SampleType, StainDate, SlideID, Sample, FOV, TotalSpores, TotalSpores_LB, TotalSpores_LBcorr, 
         log_volume_offset_m3, MedianPM2.5_ug.m3, MedianPM10_ug.m3, MedianMR, 
         SmokeLevel, MedianMCE 
         )

fasmee23_spores <- fasmee23_spores %>%
  rename(SlideID = "SampleID") %>%
  select(Project, SampleType, StainDate, SlideID, Sample, FOV, TotalSpores, TotalSpores_LB, TotalSpores_LBcorr, 
         log_volume_offset_m3, MedianPM2.5_ug.m3, MedianPM10_ug.m3, MedianMR, 
         SmokeLevel, MedianMCE 
  )

fasmee23_k4_spores <- bind_rows(k4_spores, fasmee23_spores)

fasmee23_k4_spores <- fasmee23_k4_spores %>%
  mutate(
    Project = as.factor(Project),
    SampleType = as.factor(SampleType),
    Sample = as.factor(Sample),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")),
    StainDate = as.factor(StainDate),
    log1TotalSpores_LB = log1p(TotalSpores_LB)
  ) 

ggplot(fasmee23_k4_spores, aes(x = SmokeLevel, y = MedianMR)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw()


# TotalSpores_m <- glmmTMB(TotalSpores_LBcorr ~ SmokeLevel + MedianMR + Project + offset(log_volume_offset_m3) + 
#                            (1|Sample),
#                          family=tweedie(link="log"), dispformula = ~SmokeLevel, data = fasmee23_k4_spores, ziformula = ~0)
# summary(TotalSpores_m)


TotalSpores_m <- glmmTMB(TotalSpores ~ SmokeLevel + Project + MedianMR + offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) +
                           (1|Sample),
                         family=nbinom2(link="log"), dispformula = ~SmokeLevel, data = fasmee23_k4_spores, ziformula = ~SampleType)
summary(TotalSpores_m)

simulationOutput <- simulateResiduals(fittedModel = TotalSpores_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

# Smoke spores
#-------------------------------------------------------------------------------------------------------------------------------------

smoke_fasmee23_k4_spores <- fasmee23_k4_spores %>% 
  filter(SampleType == "Smoke") 

smoke_fasmee23_k4_spores$logPM25 <- log(smoke_fasmee23_k4_spores$MedianPM2.5_ug.m3)
smoke_fasmee23_k4_spores$logPM10 <- log(smoke_fasmee23_k4_spores$MedianPM10_ug.m3)


# model_logPM <- glmmTMB(TotalSpores_LBcorr ~ logPM25*MedianMCE + Project + MedianMR +
#                          offset(log_volume_offset_m3) + (1|Sample),
#                        family=tweedie(link="log"), data = smoke_fasmee23_k4_spores)
# summary(model_logPM)

# model_logPM <- glmmTMB(TotalSpores_LBcorr ~ logPM25 + MedianMCE + MedianMR + Project +
#                          logPM25:MedianMCE + logPM25:MedianMR + Project:MedianMCE + 
#                          offset(log_volume_offset_m3) + (1|Sample),
#                        family=tweedie(link="log"), ziformula = ~0, data = smoke_fasmee23_k4_spores)
# 
# summary(model_logPM)

model_logPM <- glmmTMB(TotalSpores ~ logPM25*MedianMCE*Project + 
                         offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) + (1|Sample),
                       family=nbinom2(link="log"), dispformula = ~logPM25*MedianMCE, ziformula = ~0, data = smoke_fasmee23_k4_spores)

summary(model_logPM)

model_logPM <- glmmTMB(TotalSpores ~ SmokeLevel*Project + 
                         offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) + (1|Sample),
                       family=nbinom2(link="log"), ziformula = ~0, data = smoke_fasmee23_k4_spores)

summary(model_logPM)

simulationOutput <- simulateResiduals(fittedModel = model_logPM, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)



