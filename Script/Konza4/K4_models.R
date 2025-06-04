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
library(tweedie)
library(splines)
# library(mgcv)
# library(vegan)
#library(ggstatsplot)

spores_pa_C <- read_csv('./Output/Output_data/K4/k4_Spores_PA_C.csv') %>%
  mutate(
    Project = as.factor(Project),
    SampleType = as.factor(SampleType),
    Day = as.factor(Day),
    LB_Batch = as.factor(LB_Batch))

smoke_spores <- spores_pa_C %>%
  filter(SampleType == "Smoke")

#Konza 4 Fungi ##################################################################################################################

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


power <- 1.1265
psi <- qlogis(power - 1.0)

SmokeAmbient_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SampleType + (1|SampleID) + (1|LB_Batch),
                         family=tweedie(link="log"), data = spores_pa_C, 
                         ziformula = ~-1, start = list(psi = psi), map = list(psi = factor(NA)),
                         dispformula = ~SampleType)
summary(SmokeAmbient_model)

family_params(SmokeAmbient_model)

simulationOutput <- simulateResiduals(fittedModel = SmokeAmbient_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)


emm_log <- emmeans(SmokeAmbient_model, ~ SampleType)
emm_log

emm_response <- emmeans(SmokeAmbient_model, ~ SampleType, type = "response")
emm_response

rate_ratios <- pairs(emm_log, reverse = TRUE, type = "response")
rate_ratios

#Smoke model
#---------------------------------------------------------------------------------------------------------------------------------
smoke_spores$MeanPM2.5_ug.m3

power <- 1.1
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ log(MedianPM2.5_ug.m3) + MaxTemp_C + MeanRH + (1|SampleID) + (1|LB_Batch),
                              family=tweedie(link="log"), data = smoke_spores, 
                              ziformula = ~-1, start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_model)

family_params(Smoke_model)

simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

em_spores_SmokeLevel <- emmeans(SmokeAmbient_model, ~logPM2.5, at = list(logPM2.5 = 1.3, MedianMR = 4.8), type = "response")

summary(em_spores_SmokeLevel, infer = T)


em_contrast <- emmeans(SmokeAmbient_model, ~ logPM2.5,
                       at = list(logPM2.5 = c(1.3, 7), 
                                 MedianMR = mean(spores_pa_C$MedianMR)),
                       type = "response")

contrast(em_contrast, method = "revpairwise", infer = TRUE)

medianMR_vals <- spores_pa_C %>%
  summarise(min = min(MedianMR, na.rm = TRUE),
            mean = mean(MedianMR, na.rm = TRUE),
            max = max(MedianMR, na.rm = TRUE))

em_grid <- emmeans(SmokeAmbient_model, 
                   ~ logPM2.5 | MedianMR,
                   at = list(
                     logPM2.5 = seq(min(spores_pa_C$logPM2.5), max(spores_pa_C$logPM2.5), length.out = 100),
                     MedianMR = c(medianMR_vals$min, medianMR_vals$mean, medianMR_vals$max)
                   ),
                   type = "response")

em_df <- as.data.frame(em_grid)
em_df$PM2.5 <- exp(em_df$logPM2.5)
em_df$MR_level <- factor(em_df$MedianMR, 
                         labels = c("Min MR", "Mean MR", "Max MR"))

ggplot() +
  # Raw data
  geom_point(data = spores_pa_C, 
             aes(x = exp(logPM2.5), y = TotalSpores_LBcorr_m3), 
             alpha = 0.3, color = "gray50") +
  
  # Confidence ribbons
  geom_ribbon(data = em_df, 
              aes(x = PM2.5, ymin = asymp.LCL, ymax = asymp.UCL, fill = MR_level), 
              alpha = 0.2) +
  
  # Prediction lines
  geom_line(data = em_df, 
            aes(x = PM2.5, y = response, color = MR_level), 
            size = 1.2) +
  
  labs(
    x = expression("PM"[2.5]*" ("*mu*"g/m³)"),
    y = "Predicted Total Spores (per m³)",
    color = "Median MR Level",
    fill = "Median MR Level"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme(legend.position = "top")


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

power <- 1.45
psi <- qlogis(power - 1.0)
emission_model <- glmmTMB(spores.kg ~ poly(MeanMCE, 2) + (1|SampleID),
                       family=tweedie(link="log"), data = emission_spores, 
                       ziformula = ~poly(MeanMCE, 2))
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

# Bacteria smoke model
#-------------------------------------------------------------------------------------------------------------------------

smoke_bacteria <- read_csv('./Output/Output_data/K4/k4_Bacteria_PA_C.csv') %>%
  filter(SampleType == "Smoke") %>%
  mutate(
    SampleType = as.factor(SampleType),
    Day = as.factor(Day),
    LB_Batch = as.factor(LB_Batch))

power <- 1.1
psi <- qlogis(power - 1.0)

Smoke_BacteriaModel <- glmmTMB(TotalCells_Bcorr.m3 ~ logPM2.5 + MeanTemp_C + (1|SampleID) + (1|LB_Batch),
                       family=tweedie(link="log"), data = smoke_bacteria, 
                       ziformula = ~-1) #, start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_BacteriaModel)

family_params(Smoke_BacteriaModel)

simulationOutput <- simulateResiduals(fittedModel = Smoke_BacteriaModel, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)



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









