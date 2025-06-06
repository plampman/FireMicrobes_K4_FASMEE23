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
#library(lme4)
#library(nlme)
library(glmmTMB)
library(DHARMa)
library(tweedie)
library(splines)
# library(mgcv)
# library(vegan)
#library(ggstatsplot)

#Spores
#------------------------------------------------------------------------------------------

k4_spores <- read_csv('./Output/Output_data/K4/k4_Spores_PA_C.csv') %>%
  select(TotalSpores_FBLBcorr.m3, TotalSpores_Bcorr.m3, spores.mg, biomass_mg, SampleType, Project, Day, SampleID, 
         log_volume_offset_m3, SmokeLevel, MedianPM2.5_ug.m3, MeanPM2.5_ug.m3, LB_Batch, MeanTemp_C, 
         MeanMR, MaxTemp_C, MeanRH, MedianMCE)
fasmee23_spores <- read.csv('./Output/Output_data/FASMEE23/FASMEE23_Spores_PA_C.csv') %>%
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient")) %>%
  select(TotalSpores_FBLBcorr.m3, TotalSpores_Bcorr.m3, spores.mg, biomass_mg, SampleType, Project, Day, SampleID, 
         log_volume_offset_m3, SmokeLevel, MedianPM2.5_ug.m3, MeanPM2.5_ug.m3, LB_Batch, MeanTemp_C, 
         MeanMR, MaxTemp_C, MeanRH, MedianMCE)

fasmee23_k4_spores <- bind_rows(k4_spores, fasmee23_spores)

fasmee23_k4_spores <- fasmee23_k4_spores %>%
  mutate(
    Project = as.factor(Project),
    SampleType = as.factor(SampleType),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")),
    Day = if_else(Project == "FASMEE23", paste0("F23", "_", Day), paste0("K4", "_", Day)),
    Day = as.factor(Day),
    LB_Batch = as.factor(LB_Batch)) 

smoke_spores <- fasmee23_k4_spores %>%
  filter(SampleType == "Smoke")

emission_spores <- smoke_spores %>% filter(!is.na(MedianMCE)) 

# Smoke vs. Ambient 
#---------------------------------------------------------------------------------------------------------------

tweedie_profile <- tweedie.profile(
  TotalSpores_LBcorr_m3 ~ SampleType,
  data = spores_pa_C,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max


power <- 1.20
psi <- qlogis(power - 1.0)

SmokeAmbient_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SampleType + Project + (1|SampleID) + (1|LB_Batch),
                              family=tweedie(link="log"), data = fasmee23_k4_spores, 
                              ziformula = ~-1, #start = list(psi = psi), map = list(psi = factor(NA)),
                              dispformula = ~ Project + SampleType)  
summary(SmokeAmbient_model)


SmokeAmbient_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SmokeLevel + Project + (1|SampleID) + (1|LB_Batch),
                              family=tweedie(link="log"), data = fasmee23_k4_spores, 
                              ziformula = ~-1, dispformula = ~ Project + SmokeLevel)
summary(SmokeAmbient_model)

family_params(SmokeAmbient_model)

simulationOutput <- simulateResiduals(fittedModel = SmokeAmbient_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

#Smoke
#-------------------------------------------------------------------------------------------------------

power <- 1.1
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalSpores_Bcorr.m3 ~ Project*log(MedianPM2.5_ug.m3) + (1|SampleID) + (1|LB_Batch),
                       family=tweedie(link="log"), data = smoke_spores, 
                       ziformula = ~-1,
                       dispformula = ~log(MedianPM2.5_ug.m3)) #start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_model)

family_params(Smoke_model)

simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


# Emission model
#----------------------------------------------------------------------------------------------------

ggplot(emission_spores, aes(x = MeanMCE, y = spores.kg, color = SampleID)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1)) +
  theme_minimal()


tweedie_profile <- tweedie.profile(
  spores.kg ~ TotalSpores_Bcorr.m3 ~ Project + log(MeanPM2.5_ug.m3),
  data = emission_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.9, by = 0.01)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.27
psi <- qlogis(power - 1.0)
emission_model <- glmmTMB(TotalSpores_Bcorr.m3 ~ Project*log(MedianPM2.5_ug.m3) + offset(log(biomass_mg)) + (1|SampleID) + (1|LB_Batch),
                          family=tweedie(link="log"), data = emission_spores, 
                          ziformula = ~-1) #, start = list(psi = psi), map = list(psi = factor(NA)))

summary(emission_model)

family_params(emission_model)

simulationOutput <- simulateResiduals(fittedModel = emission_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


