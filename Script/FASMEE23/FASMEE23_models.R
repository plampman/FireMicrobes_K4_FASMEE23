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

family_params(SmokeAmbient_model)

simulationOutput <- simulateResiduals(fittedModel = SmokeAmbient_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

#Smoke model
#---------------------------------------------------------------------------------------------------------------------------------

smoke_spore_stats <- smoke_spores %>%
  group_by(SampleID, LB_Batch) %>%
  summarise(
    RH = mean(MeanRH),
    PM2.5 = mean(MedianPM2.5_ug.m3),
    Temp = mean(MeanTemp_C),
    FOVs = n(),
    NAFOVs = sum(is.na(FOV)),
    Bcorr_spores = mean(TotalSpores_Bcorr),
    Bcorr_spores_med = median(TotalSpores_Bcorr),
    volume = mean(RepVolume_L),
    MCE = mean(MeanMCE)
  )

tweedie_profile <- tweedie.profile(
  TotalSpores_Bcorr.m3 ~ logPM2.5 + MaxTemp_C + MeanRH,
  data = smoke_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max


power <- 1.44
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalSpores_Bcorr.m3 ~ logPM2.5 + MeanRH + MaxTemp_C + (1|SampleID) + (1|LB_Batch),
                       family=tweedie(link="log"), data = smoke_spores, 
                       ziformula = ~-, dispformula = ~logPM2.5) #start = list(psi = psi), map = list(psi = factor(NA)))
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


# Total cell model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------

bacteria_stat_test <- bacteria_stat_test %>% filter(SampleID != "B1B")

tweedie_profile <- tweedie.profile(
  TotalCells_FBLBcorr.m3 ~ SampleType,
  data = bacteria_stat_test,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.9, by = 0.01)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

family_params(TotalCells_m)


power <- 1.69
psi <- qlogis(power - 1.0)

TotalCells_m <- glmmTMB(TotalCells_FBLBcorr.m3 ~ SampleType + (1|SampleID) + (1|LB_Batch) + (1|Platform),
                      family=tweedie(link="log"), data = bacteria_stat_test, 
                      ziformula = ~-1, dispformula = ~SampleType)
                      #start = list(psi = psi), map = list(psi = factor(NA))))
summary(TotalCells_m)


simulationOutput <- simulateResiduals(fittedModel = TotalCells_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testZeroInflation(simulationOutput)




# TotalCells_m <- glmmTMB(TotalCells_LBcorr ~ SmokeLevel*MedianMR  + offset(log_volume_offset_m3) + 
#                           (1|SampleID),
#                         family=tweedie(link="log"), data = bacteria_blue_pa, ziformula = ~SampleType)
# summary(TotalCells_m)

family_params(TotalCells_m)

optimal_power <- 1.686
# Calculate psi value (internal parameter)
psi <- qlogis(optimal_power - 1.0)
# Fit model with fixed power parameter
TotalCells_m <- glmmTMB(
  TotalCells_LBcorr ~ SmokeLevel*MedianMR + 
    offset(log_volume_offset_m3) + 
    (1|SampleID),
  family = tweedie(link = "log"),
  data = bacteria_blue_pa, 
  ziformula = ~SampleType,
  start = list(psi = psi),  # Set the starting value 
  map = list(psi = factor(NA))  # Map to fix the parameter
)

summary(TotalCells_m)

# Check the resulting power parameter
family_params(TotalCells_m)


simulationOutput <- simulateResiduals(fittedModel = TotalCells_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testZeroInflation(simulationOutput)

testDispersion(TotalCells_m)

ggplot(bacteria_blue_pa, aes(x=MedianMR, y=TotalCells_LBcorr, color=SmokeLevel)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess") +
  theme_minimal() +
  facet_wrap(~SmokeLevel)

# Create a grid of MedianMR values covering your observed range (assuming 3.0-3.8)
mr_seq <- seq(3.0, 3.8, by = 0.2)

# Get predicted values at each combination of SmokeLevel and MedianMR
em_interact <- emmeans(TotalCells_m, 
                       ~ SmokeLevel | MedianMR, 
                       at = list(MedianMR = mr_seq),
                       type = "response")

# View the results
summary(em_interact)

# Create contrasts between smoke levels at each MedianMR value
contrasts_by_mr <- contrast(em_interact, 
                            method = "trt.vs.ctrl", 
                            ref = "None", 
                            by = "MedianMR", 
                            adjust = "none")  # No adjustment within each MR value

summary(contrasts_by_mr, infer = TRUE)

# Total cell hurdle model for smoke level
#----------------------------------------------------------------------------------------------------

presence_model <- glmmTMB(presence ~ SmokeLevel*MedianMR + 
                            (1|LB_Batch:SampleID),
                          family = binomial(link = "logit"), 
                          data = bacteria_blue_pa)
summary(presence_model)

positive_only <- subset(bacteria_blue_pa, TotalCells_LBcorr > 0)

# Fit a model for positive values only
positive_model <- glmmTMB(TotalCells_LBcorr ~ SmokeLevel+MedianMR + 
                            offset(log_volume_offset_m3) + 
                            (1|LB_Batch:SampleID),
                          family = Gamma(link = "log"), 
                          data = positive_only)
summary(positive_model)

simulationOutput <- simulateResiduals(fittedModel = positive_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testDispersion(TotalCells_m)


# Live cell model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------

LiveCells_m <- glmmTMB(LiveCells_LBcorr ~ SampleType*MedianMR + offset(log_volume_offset_m3) + 
                          (1|LB_Batch:SampleID),
                        family=tweedie(link="log"), data = bacteria_blue_pa, ziformula = ~SampleType)
summary(LiveCells_m)

simulationOutput <- simulateResiduals(fittedModel = Live.Total_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testDispersion(TotalCells_m)

# Live to total cell model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------

Live.Total_m <- glmmTMB(asin(sqrt(LiveCells.TotalCells))  ~ SampleType*MedianMR + (1|LB_Batch:SampleID),
                        family=gaussian, data = bacteria_blue_pa)
summary(Live.Total_m)

simulationOutput <- simulateResiduals(fittedModel = Live.Total_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testDispersion(TotalCells_m)


# Non-linear approach
#----------------------------------------------------------------------------------------------------

model <- glmmTMB(TotalCells_LBcorr ~ poly(MedianMR, 2)*SampleType + 
                   offset(log_volume_offset_m3) + 
                   (1|LB_Batch:SampleID),
                 ziformula = ~SampleType,
                 family = tweedie(),
                 data = bacteria_blue_pa)
summary(model)

