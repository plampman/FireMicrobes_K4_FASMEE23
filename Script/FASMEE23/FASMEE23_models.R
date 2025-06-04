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
    LB_Batch = as.factor(LB_Batch)) 

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

SmokeAmbient_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SampleType + (1|SampleID) + (1|LB_Batch),
                              family=tweedie(link="log"), data = spores_pa_C, 
                              ziformula = ~-1, dispformula = ~SampleType)
summary(SmokeAmbient_model)

emm_log <- emmeans(SmokeAmbient_model, ~ SampleType)
emm_log

emm_response <- emmeans(SmokeAmbient_model, ~ SampleType, type = "response")
emm_response

rate_ratios <- pairs(emm_log, reverse = TRUE, type = "response")
rate_ratios

family_params(SmokeAmbient_model)

simulationOutput <- simulateResiduals(fittedModel = SmokeAmbient_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

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

Smoke_model <- glmmTMB(TotalSpores_Bcorr.m3 ~ poly(logPM2.5, 3) + (1|SampleID),
                       family=tweedie(link="log"), data = smoke_spores, 
                       ziformula = ~-1, dispformula = ~logPM2.5)#, start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_model)

family_params(Smoke_model)

simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)



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
                          ziformula = ~-1, 
                          dispformula = ~ poly(MeanMCE, 2),
                          start = list(psi = psi), map = list(psi = factor(NA)))
summary(emission_model)

family_params(emission_model)

simulationOutput <- simulateResiduals(fittedModel = emission_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


library(ggeffects)

# Get predictions
pred <- ggpredict(emission_model, terms = "MeanMCE [all]")

# Plot
p <- ggplot(emission_spores, aes(x = MeanMCE, y = spores.kg)) +
  geom_point(aes(color = SampleID), alpha = 0.6, size = 2) +
  geom_ribbon(data = pred, aes(x = x, y = predicted, 
                               ymin = conf.low, ymax = conf.high),
              alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = pred, aes(x = x, y = predicted), 
            size = 1.2, color = "blue", inherit.aes = FALSE) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1e0, 1e2, 1e4, 1e6, 1e8),
                     labels = scales::scientific) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "Mean MCE", 
       y = expression(paste("Spores kg"^-1)),
       title = "Non-linear relationship between MCE and spore emission")

print(p)




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
    LB_Batch = as.factor(LB_Batch))

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

simulationOutput <- simulateResiduals(fittedModel = SmokeAmbient_BacteriaModel, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testZeroInflation(simulationOutput)


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

simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

  



