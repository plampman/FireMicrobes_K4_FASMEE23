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
    SmokeLevel = if_else(Platform == "Red", "None", SmokeLevel),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "High"))) 

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

#Spore smoke category model
#--------------------------------------------------------------------------------------------------------------------------------


# SporeCatModel <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SmokeLevel + (1|SampleID) + (1|Platform),
#                          family=tweedie(link="log"), data = spores_pa_C, dispformula = ~SmokeLevel) 
# 
# summary(SporeCatModel)


SporeCatModel <- glmmTMB(TotalSpores ~ SmokeLevel + offset(log(RepVolume_L)) + 
                                (1|Platform:SampleID) + (1|LB_Batch),
                              family=nbinom2(link="log"), data = spores_pa_C, 
                              ziformula = ~1, dispformula = ~SmokeLevel)
summary(SporeCatModel)


family_params(SporeCatModel)

res <- simulateResiduals(fittedModel = SporeCatModel, plot = F)
plotResiduals(res)
plotQQunif(res)

plotResiduals(res, spores_pa_C$SmokeLevel)  # Residuals by group
testZeroInflation(res)
testDispersion(res) 


emm_log <- emmeans(SporeCatModel, ~ SmokeLevel)
emm_log

emm_response <- emmeans(SporeCatModel, ~ SmokeLevel, type = "response")
emm_response

rate_ratios_vs_none <- contrast(emm_log, method = "trt.vs.ctrl", ref = "None", type = "response")
rate_ratios_vs_none

spore_cat <- as.data.frame(emm_response) %>% 
  mutate(
    Spores.m3 = (response*7225)/0.03,
    CI_lower = (asymp.LCL * 7225) / 0.03,
    CI_upper = (asymp.UCL * 7225) / 0.03
  )


# Create contrast dataframe from your rate_ratios_vs_none output
contrast_df <- data.frame(
  contrast = c("Low / None", "High / None"),
  ratio = c(6.06, 3.02),
  SE = c(2.41, 1.27),
  p.value = c(0.0001, 0.0161),  # <.0001 converted to 0.0001
  stringsAsFactors = FALSE
)

# Create significance bars 
sig_bars <- data.frame(
  x = c(1, 1),  # Both start from position 1 (None)
  xend = c(2, 3),  # End at positions 2 (Low) and 3 (High)
  pair = c("Low / None", "High / None"),
  stringsAsFactors = FALSE
)

# Extract p-values using contrast names
sig_bars$p_value <- sapply(sig_bars$pair, function(pair) {
  idx <- which(contrast_df$contrast == pair)
  if(length(idx) > 0) {
    return(contrast_df$p.value[idx])
  } else {
    return(1.0) 
  }
})

# Calculate maximum height for positioning significance bars
max_height <- max(spore_cat$CI_upper, na.rm = TRUE)

# Position bars at different heights to avoid overlap
sig_bars$y <- max_height * c(3.8, 3.9)

# Create significance labels
sig_bars$label <- vapply(sig_bars$p_value, function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}, character(1))

# Create the plot
p1 <- ggplot(spore_cat, aes(x = SmokeLevel, y = Spores.m3)) +
  geom_col(fill = "#0072B2", alpha = 0.6) +
  # Add raw data points with jitter (filter out zeros and NAs for log scale)
  geom_point(data = spores_pa_C, 
             aes(x = SmokeLevel, y = TotalSpores_FBLBcorr.m3), 
             position = position_jitter(width = 0.1), 
             alpha = 0.5, size = 2, color = "#0072B2") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, color = "black", linewidth = 0.7) +
  labs(x = expression(paste("PM"[2.5], " μg/m"^3, " range")), 
       y = expression(paste("Spores/m"^3)),
       title = expression(paste("Utah Subalpine Fungal Spores by ", "PM"[2.5], " Range"))) +
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  # Add significance bars
  geom_segment(data = sig_bars, 
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.7) +
  geom_segment(data = sig_bars,
               aes(x = x, xend = x, y = y, yend = y - max_height * 0.02),
               linewidth = 0.7) +
  geom_segment(data = sig_bars,
               aes(x = xend, xend = xend, y = y, yend = y - max_height * 0.02),
               linewidth = 0.7) +
  geom_text(data = sig_bars,
            aes(x = (x + xend) / 2, y = y + max_height * 0.03),
            label = sig_bars$label,
            size = 3.5) +
  # Customize x-axis labels
  scale_x_discrete(labels = c("None\n(<20)", "Low\n(20-400)", "High\n(≥400)")) +
  coord_cartesian(ylim = c(0, max(spores_pa_C$TotalSpores_FBLBcorr.m3 * 0.6)))

print(p1)


ggsave("./Output/Output_figs/FASMEE23/FASMEE23_SmokeLevel_Spores.png", p1, width = 5, height = 5, dpi = 600)

#Smoke model
#---------------------------------------------------------------------------------------------------------------------------------

smoke_spores$MedianPM10_ug.m3

ggplot(smoke_spores, aes(x = MedianPM10_ug.m3, y = TotalSpores_FBLBcorr.m3, color = SampleID)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1)) +
  theme_minimal()

tweedie_profile <- tweedie.profile(
  TotalSpores.m3 ~ poly(log(MedianPM2.5_ug.m3), 2) + TotalSpores_FBLBcorr,
  data = smoke_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max


power <- 1.005647
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalSpores.m3 ~ log(MedianPM2.5_ug.m3) + TotalSpores_FBLBcorr +
                         (1|SampleID),
                       family=tweedie(link="log"), data = smoke_spores, 
                       ziformula = ~1, dispformula = ~log(MedianPM2.5_ug.m3),
                       start = list(psi = psi), map = list(psi = factor(NA)))
                       
summary(Smoke_model)

family_params(Smoke_model)


#, start = list(psi = psi), map = list(psi = factor(NA)))

unique(smoke_spores$TotalSpores_FBLB)

Smoke_model <- glmmTMB(TotalSpores ~ poly(log(MedianPM2.5_ug.m3), 3) + offset(log(RepVolume_m3)) + 
                         (1|SampleID) + (1|LB_Batch),
                       family=nbinom2(link="log"), data = smoke_spores, 
                       ziformula = ~1, dispformula = ~poly(log(MedianPM2.5_ug.m3), 3))
summary(Smoke_model)


res <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotResiduals(res)
plotQQunif(res)

testZeroInflation(res)
testDispersion(res) 


pm25_emmeans <- emmeans(Smoke_model, 
                        ~MedianPM2.5_ug.m3,  
                        at = list(MedianPM2.5_ug.m3 = seq(min(smoke_spores$MedianPM2.5_ug.m3, na.rm = TRUE),
                                                          max(smoke_spores$MedianPM2.5_ug.m3, na.rm = TRUE),
                                                          length.out = 50)),
                        type = "response")

pm25_plot_data <- as.data.frame(pm25_emmeans)

pm25_values <- pm25_plot_data$MedianPM2.5_ug.m3
predictions <- pm25_plot_data$response

cubic_fit <- lm(predictions ~ poly(log(pm25_values), 3, raw = TRUE))
coeffs <- coef(cubic_fit)

# Your calculator formula
cat("TotalSpores =", coeffs[1], "+", coeffs[2], "* log(PM2.5) +", 
    coeffs[3], "* [log(PM2.5)]^2 +", coeffs[4], "* [log(PM2.5)]^3")


p2 <- ggplot() +
  geom_point(data = smoke_spores,
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores.m3),
             alpha = 0.6, size = 2, color = "black") +
  geom_ribbon(data = pm25_plot_data, 
              aes(x = MedianPM2.5_ug.m3, ymin = (asymp.LCL*7225)/0.03, ymax = (asymp.UCL*7225)/0.03), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = pm25_plot_data, 
            aes(x = MedianPM2.5_ug.m3, y = (response*7225)/0.03), 
            size = 1.2, color = "blue") +
  labs(x = expression(PM[2.5]~(μg/m^3)),
       y = expression(Spores~m^-3),
       title = expression(paste("Utah Subalpine Smoke Spore concentration vs ", "PM"[2.5]))) +
  coord_cartesian(ylim = c(0, max(spores_pa_C$TotalSpores_FBLBcorr.m3 * 0.6))) +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0, face = "bold")) +
  theme_minimal()


ggsave("./Output/Output_figs/FASMEE23/FASMEE23_Spores_PM2.5.png", p2, width = 5.5, height = 5, dpi = 600, bg = "white")


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
    SmokeLevel = if_else(Platform == "Red", "None", SmokeLevel),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "High"))) 

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

#Bacteria smoke category model
#--------------------------------------------------------------------------------------------------------------------------------


BacteriaCatModel <- glmmTMB(TotalCells ~ SmokeLevel + offset(log(RepVolume_L)) + 
                           (1|Platform:SampleID) + (1|LB_Batch),
                         family=nbinom2(link="log"), data = bacteria, 
                         ziformula = ~1, dispformula = ~SmokeLevel)
summary(BacteriaCatModel)


res <- simulateResiduals(fittedModel = BacteriaCatModel, plot = F)
plotResiduals(res)
plotQQunif(res)

plotResiduals(res, bacteria$SmokeLevel)  # Residuals by group
testZeroInflation(res)
testDispersion(res) 


emm_log <- emmeans(BacteriaCatModel, ~ SmokeLevel)
emm_log

emm_response <- emmeans(BacteriaCatModel, ~ SmokeLevel, type = "response")
emm_response

rate_ratios_vs_none <- contrast(emm_log, method = "trt.vs.ctrl", ref = "None", type = "response")
rate_ratios_vs_none

bacteria_cat <- as.data.frame(emm_response) %>% 
  mutate(
    Bacteria.m3 = (response*7225)/0.03,
    CI_lower = (asymp.LCL * 7225) / 0.03,
    CI_upper = (asymp.UCL * 7225) / 0.03
  )


# Create contrast dataframe from your rate_ratios_vs_none output
contrast_df <- data.frame(
  contrast = c("Low / None", "High / None"),
  ratio = c(6.06, 3.02),
  SE = c(2.41, 1.27),
  p.value = c(0.049, 0.3049),  # <.0001 converted to 0.0001
  stringsAsFactors = FALSE
)

# Create significance bars 
sig_bars <- data.frame(
  x = c(1, 1),  # Both start from position 1 (None)
  xend = c(2, 3),  # End at positions 2 (Low) and 3 (High)
  pair = c("Low / None", "High / None"),
  stringsAsFactors = FALSE
)

# Extract p-values using contrast names
sig_bars$p_value <- sapply(sig_bars$pair, function(pair) {
  idx <- which(contrast_df$contrast == pair)
  if(length(idx) > 0) {
    return(contrast_df$p.value[idx])
  } else {
    return(1.0) 
  }
})

# Calculate maximum height for positioning significance bars
max_height <- max(bacteria_cat$CI_upper, na.rm = TRUE)

# Position bars at different heights to avoid overlap
sig_bars$y <- max_height * c(3.8, 3.95)

# Create significance labels
sig_bars$label <- vapply(sig_bars$p_value, function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}, character(1))

# Create the plot
p1 <- ggplot(bacteria_cat, aes(x = SmokeLevel, y = Bacteria.m3)) +
  geom_col(fill = "#0072B2", alpha = 0.6) +
  # Add raw data points with jitter (filter out zeros and NAs for log scale)
  geom_point(data = bacteria, 
             aes(x = SmokeLevel, y = TotalCells_FBLBcorr.m3), 
             position = position_jitter(width = 0.1), 
             alpha = 0.5, size = 2, color = "#0072B2") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, color = "black", linewidth = 0.7) +
  labs(x = expression(paste("PM"[2.5], " μg/m"^3, " range")), 
       y = expression(paste("Bacteria/m"^3)),
       title = expression(paste("Utah Subalpine Bacteria by ", "PM"[2.5], " Range"))) +
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  # Add significance bars
  geom_segment(data = sig_bars, 
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.7) +
  geom_segment(data = sig_bars,
               aes(x = x, xend = x, y = y, yend = y - max_height * 0.03),
               linewidth = 0.7) +
  geom_segment(data = sig_bars,
               aes(x = xend, xend = xend, y = y, yend = y - max_height * 0.03),
               linewidth = 0.7) +
  geom_text(data = sig_bars,
            aes(x = (x + xend) / 2, y = y + max_height * 0.07),
            label = sig_bars$label,
            size = 3.5) +
  # Customize x-axis labels
  scale_x_discrete(labels = c("None\n(<20)", "Low\n(20-400)", "High\n(≥400)")) +
  coord_cartesian(ylim = c(0, max(bacteria$TotalCells_FBLBcorr.m3 * 0.35)))

print(p1)


ggsave("./Output/Output_figs/FASMEE23/FASMEE23_SmokeLevel_Bacteria.png", p1, width = 5, height = 5, dpi = 600)


# Smoke only model
#----------------------------------------------------------------------------------------------------------------------------------

ggplot(smoke_bacteria, aes(x = logPM2.5, y = TotalCells_Bcorr.m3, color = SampleID)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1)) +
  theme_minimal()


tweedie_profile <- tweedie.profile(
  TotalCells.m3 ~ poly(log(MedianPM2.5_ug.m3), 2) + TotalCells_FBLBcorr, 
  data = smoke_bacteria,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.4
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalCells.m3 ~ poly(log(MedianPM2.5_ug.m3), 2) + TotalCells_FBLBcorr +
                         (1|SampleID),
                       family=tweedie(link="log"), data = smoke_bacteria, 
                       ziformula = ~1, dispformula = ~logPM2.5,
                       start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_model)

family_params(Smoke_model)

Smoke_model <- glmmTMB(TotalCells ~ poly(log(MedianPM2.5_ug.m3), 4) + TotalCells_FBLB + offset(log(RepVolume_m3)) + 
                         (1|SampleID),
                       family=nbinom2(link="log"), data = smoke_bacteria, 
                       ziformula = ~1, dispformula = ~poly(log(MedianPM2.5_ug.m3), 4))
summary(Smoke_model)

simulationOutput <- simulateResiduals(fittedModel = Smoke_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

  
pm25_emmeans <- emmeans(Smoke_model, 
                        ~MedianPM2.5_ug.m3,  
                        at = list(MedianPM2.5_ug.m3 = seq(min(smoke_bacteria$MedianPM2.5_ug.m3, na.rm = TRUE),
                                                          max(smoke_bacteria$MedianPM2.5_ug.m3, na.rm = TRUE),
                                                          length.out = 50),
                                  TotalSpores_FBLBcorr = 0),
                        type = "response")

pm25_plot_data <- as.data.frame(pm25_emmeans)

pm25_values <- pm25_plot_data$MedianPM2.5_ug.m3
predictions <- pm25_plot_data$response

cubic_fit <- lm(predictions ~ poly(log(pm25_values), 4, raw = TRUE))
coeffs <- coef(cubic_fit)

# Your calculator formula
cat("Bacteria/m3 =", coeffs[1], "+", coeffs[2], "* log(PM2.5) +", 
    coeffs[3], "* [log(PM2.5)]^2 +", coeffs[4], "* [log(PM2.5)]^3 +",
    coeffs[5], "* [log(PM2.5)]^4")


p2 <- ggplot() +
  geom_point(data = smoke_bacteria,
             aes(x = MedianPM2.5_ug.m3, y = TotalCells.m3),
             alpha = 0.6, size = 2, color = "black") +
  geom_ribbon(data = pm25_plot_data, 
              aes(x = MedianPM2.5_ug.m3, ymin = (asymp.LCL*7225)/0.03, ymax = (asymp.UCL*7225)/0.03), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = pm25_plot_data, 
            aes(x = MedianPM2.5_ug.m3, y = (response*7225)/0.03), 
            size = 1.2, color = "blue") +
  labs(x = expression(PM[2.5]~(μg/m^3)),
       y = expression(Bacteria~m^-3),
       title = expression(paste("Utah Subalpine Smoke Bacteria concentration vs ", "PM"[2.5]))) +
  coord_cartesian(ylim = c(0, max(smoke_bacteria$TotalCells.m3 * 0.4))) +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0, face = "bold")) +
  theme_minimal()


ggsave("./Output/Output_figs/FASMEE23/FASMEE23_Bacteria_PM2.5.png", p2, width = 5.5, height = 5, dpi = 600, bg = "white")



