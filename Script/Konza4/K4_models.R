#########################################################################
###
### Program name: Konza 4 statistics
###
### Purpose: Konza 4 mixed models
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/15/2025
###
#########################################################################

### Run Konza4_Metadata.R, K4_PA.R, K4_Carbon.R, and K4_SporeConcentration.R first


library(tidyverse)
library(emmeans)
#library(lme4)
#library(nlme)
library(glmmTMB)
library(DHARMa)
# library(mgcv)
# library(vegan)
#library(ggstatsplot)

#Konza 4 Fungi
#--------------------------------------------------------------------------------------


TotalSpores_SmokeLevel <- glmmTMB(TotalSpores_LBcorr_m3 ~ SmokeLevel + (1|SampleID) + (1|RepVolume_m3),
                         family=tweedie(link="log"), data = spores_pa_C, ziformula = ~0)
summary(TotalSpores_SmokeLevel)

simulationOutput <- simulateResiduals(fittedModel = TotalSpores_SmokeLevel, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

em_spores_SmokeLevel <- emmeans(TotalSpores_SmokeLevel, ~ SmokeLevel, type = "response")
# View the predicted means
summary(em_spores_SmokeLevel, infer = T)

# Convert to ratios instead of differences
contrasts_ratio_spores_SmokeLevel <- contrast(em_spores_SmokeLevel, method = "trt.vs.ctrl", ref = "None", type = "response")
summary(contrasts_ratio_spores_SmokeLevel, infer = TRUE)

#No adjustment for multiple comparisons with p-values
contrasts_ratio_noadj <- contrast(em, method = "trt.vs.ctrl", ref = "None", type = "response")
summary(contrasts_ratio_noadj, infer = TRUE, adjust = "none")


# Total cell model for smoke with MCE
#----------------------------------------------------------------------------------------------------

smoke_spores_pa_C <- spores_pa_C %>% filter(SampleType == "Smoke" & !is.na(MeanMCE)) 


smoke_spores_pa_C$logPM25 <- log(smoke_spores_pa_C$MedianPM2.5_ug.m3)
smoke_spores_pa_C$logPM10 <- log(smoke_spores_pa_C$MedianPM10_ug.m3)
smoke_spores_pa_C$logPM1<- log(smoke_spores_pa_C$MedianPM1.0_ug.m3)


model_logPM <- glmmTMB(TotalSpores_LBcorr ~ MedianMCE + logPM25 + MedianMR + 
                         logPM25:MedianMCE +
                         offset(log_volume_offset_m3) + (1|Sample),
                       family=tweedie(link="log"), data = smoke_spores_pa_C)

summary(model_logPM)

slim_UI_EPA_C$MedianCO2_CO_mg.m3

model_logPM <- glmmTMB(TotalSpores_LBcorr ~ MedianMCE*logPM10 + log(MedianCO2_CO_mg.m3) +
                         offset(log_volume_offset_m3) + (1|Sample),
                       family=tweedie(link="log"), data = smoke_spores_pa_C)

summary(model_logPM)

simulationOutput <- simulateResiduals(fittedModel = model_logPM, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


# Total cell hurdle model for smoke level
#----------------------------------------------------------------------------------------------------

presence_model <- glmmTMB(presence ~ SmokeLevel + MedianMR + 
                            (1|Sample),
                          family = binomial(link = "logit"), 
                          data = spores_pa)
summary(presence_model)

positive_only <- subset(spores_pa, TotalSpores_LBcorr > 0)

# Fit a model for positive values only
positive_model <- glmmTMB(TotalSpores_LBcorr ~ SmokeLevel + MedianMR + 
                            offset(log_volume_offset_m3) + 
                            (1|Sample),
                          family = Gamma(link = "log"), 
                          data = positive_only)
summary(positive_model)

simulationOutput <- simulateResiduals(fittedModel = positive_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testDispersion(TotalCells_m)






#Extra
#--------------------------------------------------------------------------------------

#write.csv(spores_pa, './spores_pa20250406.csv', row.names = F)

spore_glm <- glm(Tot_spores.m3 ~ SampleType + MedianMR,
                         family = Gamma(link = "log"), data = spores_pa)

summary(spore_glm)

simulationOutput <- simulateResiduals(fittedModel = spore_glm, plot = F)

plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testDispersion(spore_glm)

emm <- emmeans(spore_glm, ~ SampleType, type = "response")

# View the results
emm

pairs(emmeans(spore_glm, ~ SampleType, type = "response"), reverse = TRUE, infer = TRUE)


spores_pm <- spores_pm %>%
  filter(if_all(everything(), ~ !is.na(.)))

TotalSpores_m <- glmmTMB(TotalSpores ~ SmokeLevel + 
                           offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) + (1|SampleID),
                         family=nbinom2(link="log"), data = spores_pa_C, ziformula = ~0)
summary(TotalSpores_m)


