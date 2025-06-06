#########################################################################
###
### Program name: FASMEE23 statistics
###
### Purpose: FASMEE23 mixed models
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 04/06/2025
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


spores_pa_C <- read.csv('./Output/Output_data/FASMEE23/FASMEE23_Spores_PA_C.csv') %>%
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient")) %>%
  mutate(
    Project = as.factor(Project),
    SampleType = as.factor(SampleType),
    Day = as.factor(Day),
    LB_Batch = as.factor(LB_Batch),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High"))) 

summary(spores_pa_C)

### FASMEE Spores

smoke_spores <- spores_pa_C %>%
  filter(Platform == "Blue" & SampleType == "Smoke")

# Total spore model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------

# Smoke vs. Ambient model
#---------------------------------------------------------------------------------------------------------------------------------

tweedie_profile <- tweedie.profile(
  TotalSpores_FBLBcorr.m3 ~ SampleType,
  data = spores_pa_C,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max


power <- 1.33
psi <- qlogis(power - 1.0)

SmokeAmbient_model <- glmmTMB(TotalSpores.m3 ~ SampleType + (1|SampleID) + (1|LB_Batch),
                              family=tweedie(link="log"), data = spores_pa_C, 
                              ziformula = ~-1, dispformula = ~SampleType)
summary(SmokeAmbient_model)


SmokeAmbient_model <- glmmTMB(TotalSpores ~ SampleType + offset(log(RepVolume_L)) + 
                                (1|SampleID) + (1|LB_Batch),
                              family=nbinom2(link="log"), data = spores_pa_C, 
                              ziformula = ~-1, dispformula = ~SampleType)
summary(SmokeAmbient_model)


emm_log <- emmeans(SmokeAmbient_model, ~ SampleType)
emm_log

emm_response <- emmeans(SmokeAmbient_model, ~ SampleType, type = "response")
emm_response

rate_ratios <- pairs(emm_log, reverse = TRUE, type = "response")
rate_ratios

family_params(SmokeAmbient_model)

res <- simulateResiduals(fittedModel = SmokeAmbient_model, plot = F)
plotResiduals(res)
plotQQunif(res)

plotResiduals(res, spores_pa_C$SampleType)  # Residuals by group
testZeroInflation(res)
testDispersion(res) 

#Smoke model
#---------------------------------------------------------------------------------------------------------------------------------

smoke_spores$MedianPM10_ug.m3

ggplot(smoke_spores, aes(x = MedianPM10_ug.m3, y = TotalSpores_Bcorr.m3, color = SampleID)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1)) +
  theme_minimal()

tweedie_profile <- tweedie.profile(
  TotalSpores_Bcorr.m3 ~ logPM2.5, 
  data = smoke_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.4
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalSpores.m3 ~ poly(logPM2.5, 3) + MeanTemp_C + (1|SampleID) + (1|LB_Batch),
                       family=tweedie(link="log"), data = smoke_spores, 
                       ziformula = ~-1, dispformula = ~poly(logPM2.5, 3))
summary(Smoke_model)

family_params(Smoke_model)

#, start = list(psi = psi), map = list(psi = factor(NA)))



Smoke_model <- glmmTMB(TotalSpores ~ poly(logPM2.5, 3) + MeanTemp_C + offset(log(RepVolume_m3)) + 
                         (1|SampleID) + (1|LB_Batch),
                       family=nbinom2(link="log"), data = smoke_spores, 
                       ziformula = ~-1, dispformula = ~poly(logPM2.5, 3))
summary(Smoke_model)


simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


pm25_emmeans <- emmeans(Smoke_model, 
                        ~logPM2.5,  
                        at = list(logPM2.5 = seq(min(smoke_spores$logPM2.5, na.rm = TRUE),
                                                 max(smoke_spores$logPM2.5, na.rm = TRUE),
                                                 length.out = 50)),
                        type = "response")

pm25_plot_data <- as.data.frame(pm25_emmeans)

pm25_plot_data$PM2.5 <- exp(pm25_plot_data$logPM2.5)

ggplot() +
  geom_ribbon(data = pm25_plot_data, 
              aes(x = PM2.5, ymin = (asymp.LCL*7225)/0.03, ymax = (asymp.UCL*7225)/0.03), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = pm25_plot_data, 
            aes(x = PM2.5, y = (response*7225)/0.03), 
            size = 1.2, color = "blue") +
  geom_point(data = smoke_spores, 
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores.m3), 
             alpha = 0.6, size = 2, color = "black") +
  labs(x = expression(PM[2.5]~(μg/m^3)),
       y = expression(Spores~m^-3),
       title = "Spore Concentration vs PM2.5 (Cubic Relationship)") +
  theme_minimal()


temp_emmeans <- emmeans(Smoke_model,
                        specs = "MeanTemp_C",
                        at = list(MeanTemp_C = seq(min(smoke_spores$MeanTemp_C, na.rm = TRUE),
                                                   max(smoke_spores$MeanTemp_C, na.rm = TRUE),
                                                   length.out = 30)), type = "response")

temp_plot_data <- as.data.frame(temp_emmeans)

ggplot(temp_plot_data, aes(x = MeanTemp_C, y = response)) +
  geom_line(size = 1.2, color = "red") +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
              alpha = 0.3, fill = "red") +
  labs(x = "Mean Temperature (°C)",
       y = "Predicted Total Spores",
       title = "Spore Count vs Temperature") +
  theme_minimal()



# Emission model
#----------------------------------------------------------------------------------------------------

emission_spores <- smoke_spores %>% filter(!is.na(MeanMCE)) 

ggplot(emission_spores, aes(x = MeanMCE, y = spores.kg, color = SampleID)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1)) +
  theme_minimal()


tweedie_profile <- tweedie.profile(
  spores.kg ~ poly(MeanMCE, 2),
  data = emission_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.9, by = 0.01)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.69
psi <- qlogis(power - 1.0)
emission_model <- glmmTMB(spores.kg ~ poly(MeanMCE, 2) + (1|SampleID),
                          family=tweedie(link="log"), data = emission_spores, 
                          ziformula = ~1, 
                          dispformula = ~ poly(MeanMCE, 2),
                          start = list(psi = psi), map = list(psi = factor(NA)))
summary(emission_model)

family_params(emission_model)


emission_model <- glmmTMB(TotalSpores ~ poly(MeanMCE, 2) + offset(log(RepVolume_m3)) + 
                            offset(log(biomass_Mg)) + offset(log(Ambient_Offset)) +
                            (1|SampleID) + (1|LB_Batch),
                          family=nbinom2(link="log"), data = emission_spores, 
                          ziformula = ~1)
summary(emission_model)

simulationOutput <- simulateResiduals(fittedModel = emission_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


MCE_emmeans <- emmeans(emission_model, 
                        specs = ~MeanMCE,
                        at = list(MeanMCE = seq(min(emission_spores$MeanMCE, na.rm = TRUE),
                                                 max(emission_spores$MeanMCE, na.rm = TRUE),
                                                 length.out = 50)),
                        type = "response")

MCE_plot_data <- as.data.frame(MCE_emmeans)


ggplot() +
  geom_ribbon(data = MCE_plot_data, 
              aes(x = MeanMCE, ymin = (asymp.LCL*7225)/0.03, ymax = (asymp.UCL*7225)/0.03), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = MCE_plot_data, 
            aes(x = MeanMCE, y = (response*7225)/0.03), 
            size = 1.2, color = "blue") +
  geom_point(data = emission_spores, 
             aes(x = MeanMCE, y = TotalSpores.m3), 
             alpha = 0.6, size = 2, color = "black") +
  labs(x = 'MCE',
       y = expression(Spores~m^-3),
       title = "Spore Concentration vs MCE (Cubic Relationship)") +
  theme_minimal()




########### FASMEE Bacteria ###########################################################
#######################################################################################


# Total cell model for ambient versus smoke 
#----------------------------------------------------------------------------------------------------

bacteria <- read_csv('./Output/Output_data/FASMEE23/FASMEE23_Bacteria_PA_C.csv') %>%
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient")) %>%
  mutate(
    Project = as.factor(Project),
    SampleType = as.factor(SampleType),
    Day = as.factor(Day),
    LB_Batch = as.factor(LB_Batch),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")))

smoke_bacteria <- bacteria %>%
  filter(Platform == "Blue" & SampleType == "Smoke")

tweedie_profile <- tweedie.profile(
  TotalCells_FBLBcorr.m3 ~ SampleType,
  data = bacteria,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.9, by = 0.01)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.6
psi <- qlogis(power - 1.0)

SmokeAmbient_BacteriaModel <- glmmTMB(TotalCells_FBLBcorr.m3 ~ SampleType + Platform + (1|SampleID) + (1|LB_Batch),
                      family=tweedie(link="log"), data = bacteria, 
                      ziformula = ~-1, dispformula = ~SampleType) # start = list(psi = psi), map = list(psi = factor(NA)))
summary(SmokeAmbient_BacteriaModel)

family_params(SmokeAmbient_BacteriaModel)


SmokeAmbient_BacteriaModel <- glmmTMB(TotalCells ~ SampleType + Platform + TotalCells_FBLB + offset(log(RepVolume_m3)) +
                                        (1|SampleID),
                                      family=nbinom2(link="log"), data = bacteria, 
                                      ziformula = ~1, dispformula = ~1) 
summary(SmokeAmbient_BacteriaModel)


SmokeAmbient_BacteriaModel <- glmmTMB(LIVECounts ~ SampleType + Platform + LiveCells_FBLB + offset(log(RepVolume_m3)) +
                                        (1|SampleID),
                                      family=nbinom2(link="log"), data = bacteria, 
                                      ziformula = ~1, dispformula = ~SampleType) 
summary(SmokeAmbient_BacteriaModel)

res <- simulateResiduals(fittedModel = SmokeAmbient_BacteriaModel, plot = F)
plotResiduals(res)
plotQQunif(res)

plotResiduals(res, bacteria$SampleType)  # Residuals by group
testZeroInflation(res)
testDispersion(res)


emm_log <- emmeans(SmokeAmbient_BacteriaModel, ~ SampleType)
emm_log

emm_response <- emmeans(SmokeAmbient_BacteriaModel, ~ SampleType, type = "response")
emm_response

rate_ratios <- pairs(emm_log, reverse = TRUE, type = "response")
rate_ratios

# Smoke only model
#----------------------------------------------------------------------------------------------------------------------------------

ggplot(smoke_bacteria, aes(x = logPM2.5, y = TotalCells_Bcorr.m3, color = SampleID)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1)) +
  theme_minimal()


tweedie_profile <- tweedie.profile(
  TotalSpores_Bcorr.m3 ~ logPM2.5, 
  data = smoke_bacteria,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.4
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalCells_Bcorr.m3 ~ poly(logPM2.5, 2) + MeanTemp_C + MeanRH + (1|SampleID),
                       family=tweedie(link="log"), data = smoke_bacteria, 
                       ziformula = ~-1, dispformula = ~logPM2.5)#, start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_model)

family_params(Smoke_model)

Smoke_model <- glmmTMB(TotalCells ~ poly(logPM2.5, 3) + MeanTemp_C + TotalCells_FBLB + offset(log(RepVolume_m3)) + 
                         (1|SampleID),
                       family=nbinom2(link="log"), data = smoke_bacteria, 
                       ziformula = ~1, dispformula = ~poly(logPM2.5, 3))
summary(Smoke_model)

simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

  



