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
    LB_Batch = as.factor(LB_Batch),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "High")))

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

SmokeAmbient_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SampleType + (1|SampleID),
                         family=tweedie(link="log"), data = spores_pa_C, 
                         ziformula = ~-1, start = list(psi = psi), map = list(psi = factor(NA)),
                         dispformula = ~SampleType)
summary(SmokeAmbient_model)

family_params(SmokeAmbient_model)


SmokeAmbient_model <- glmmTMB(TotalSpores ~ SampleType + TotalSpores_FBLB + offset(log(RepVolume_m3)) + 
                                (1|SampleID),
                              family=nbinom2(link="log"), data = spores_pa_C, 
                              ziformula = ~1)
summary(SmokeAmbient_model)

res <- simulateResiduals(fittedModel = SmokeAmbient_model, plot = F)
plotResiduals(res)
plotQQunif(res)

plotResiduals(res, spores_pa_C$SampleType)  # Residuals by group
testZeroInflation(res)
testDispersion(res) 


emm_log <- emmeans(SmokeAmbient_model, ~ SampleType)
emm_log

emm_response <- emmeans(SmokeAmbient_model, ~ SampleType, type = "response")
emm_response

rate_ratios <- pairs(emm_log, reverse = TRUE, type = "response")
rate_ratios

#Spore smoke category model
#--------------------------------------------------------------------------------------------------------------------------------

# 
# SporeCatModel <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ SmokeLevel + (1|SampleID),
#                               family=tweedie(link="log"), data = spores_pa_C) 

family_params(SporeCatModel)

SporeCatModel <- glmmTMB(TotalSpores ~ SmokeLevel + offset(log(Slide_RepVolume_L)) + 
                           (1|SampleID) + (1|LB_Batch),
                         family=nbinom2(link="log"), data = spores_pa_C, 
                         ziformula = ~1)
summary(SporeCatModel)


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
  p.value = c(0.0372, 0.0300),  # <.0001 converted to 0.0001
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
sig_bars$y <- max_height * c(3.1, 3.2)

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
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, color = "black", linewidth = 0.7) +
  # Add raw data points with jitter (filter out zeros and NAs for log scale)
  geom_point(data = spores_pa_C, 
             aes(x = SmokeLevel, y = TotalSpores_FBLBcorr.m3), 
             position = position_jitter(width = 0.1), 
             alpha = 0.5, size = 2, color = "#0072B2") +
  labs(x = expression(paste("PM"[2.5], " μg/m"^3, " range")), 
       y = expression(paste("Spores/m"^3)),
       title = expression(paste("Tallgrass Prairie Fungal Spores by ", "PM"[2.5], " Range"))) +
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
  coord_cartesian(ylim = c(0, max(spores_pa_C$TotalSpores_FBLBcorr.m3 * 0.5)))

print(p1)


ggsave("./Output/Output_figs/K4/Konza_SmokeLevel_Spores.png", p1, width = 5, height = 5, dpi = 600)

#Smoke model
#---------------------------------------------------------------------------------------------------------------------------------
ggplot(spores_pa_C, aes(x = MedianPM2.5_ug.m3, y = TotalSpores_FBLBcorr.m3)) +
  geom_point() +
  geom_point(data = spores_pa_C, 
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores_FBLBcorr.m3), 
             position = position_jitter(width = 0.2), 
             alpha = 0.5, size = 2, color = "black", inherit.aes = FALSE) +
  geom_smooth(method = "loess", se = FALSE, aes(group = 1), color = "blue") +
  theme_minimal()


smoke_spores <- smoke_spores %>%
  filter(SampleID != "14A")

tweedie_profile <- tweedie.profile(
  TotalSpores.m3 ~ log(MedianPM2.5_ug.m3) + TotalSpores_FBLB,
  data = smoke_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.8, by = 0.05)  # typical range: 1 < p < 2
)

tweedie_profile$p.max


power <- 1.1
psi <- qlogis(power - 1.0)

Smoke_model <- glmmTMB(TotalSpores_FBLBcorr.m3 ~ log(MedianPM2.5_ug.m3) + MaxTemp_C + MeanRH + (1|SampleID) + (1|LB_Batch),
                              family=tweedie(link="log"), data = smoke_spores, 
                       ziformula = ~-1, start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_model)

family_params(Smoke_model)


Smoke_model <- glmmTMB(TotalSpores ~ poly(log(MedianPM2.5_ug.m3), 3) + offset(log(RepVolume_m3)) + 
                         (1|SampleID) + (1|LB_Batch),
                       family=nbinom2(link="log"), data = smoke_spores, 
                       ziformula = ~1)
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

cubic_fit <- lm(predictions ~ poly(log(pm25_values), 2, raw = TRUE))
coeffs <- coef(cubic_fit)

# Your calculator formula
cat("Spores/m3 =", coeffs[1], "+", coeffs[2], "* log(PM2.5) +", 
    coeffs[3], "* [log(PM2.5)]^2")

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
       title = expression(paste("Tallgrass Prairie Smoke Spore concentration vs ", "PM"[2.5]))) +
  coord_cartesian(ylim = c(0, max(spores_pa_C$TotalSpores_FBLBcorr.m3 * 0.7))) +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0, face = "bold")) +
  theme_minimal()


ggsave("./Output/Output_figs/K4/K4_Spores_PM2.5.png", p2, width = 5.5, height = 5, dpi = 600, bg = "white")





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
  filter(SampleType == "Smoke" & SampleID != "14A") %>%
  mutate(
    SampleType = as.factor(SampleType),
    Day = as.factor(Day),
    LB_Batch = as.factor(LB_Batch)) 


tweedie_profile <- tweedie.profile(
  spores.kg ~ poly(MeanMCE, 2),
  data = emission_spores,
  method = "series",     
  do.plot = TRUE,        
  p.vec = seq(1.1, 1.9, by = 0.01)  # typical range: 1 < p < 2
)

tweedie_profile$p.max

power <- 1.1
psi <- qlogis(power - 1.0)

Smoke_BacteriaModel <- glmmTMB(TotalCells_Bcorr.m3 ~ log(MedianPM2.5_ug.m3) + (1|LB_Batch) + (1|SampleID),
                       family=tweedie(link="log"), data = smoke_bacteria, 
                       ziformula = ~1) #, start = list(psi = psi), map = list(psi = factor(NA)))
summary(Smoke_BacteriaModel)

family_params(Smoke_BacteriaModel)

Smoke_BacteriaModel <- glmmTMB(TotalCells ~ log(MedianPM2.5_ug.m3) + offset(log(RepVolume_m3)) + 
                         (1|SampleID),
                       family=nbinom2(link="log"), data = smoke_bacteria, 
                       ziformula = ~1)
summary(Smoke_BacteriaModel)

# The relationship between log(PM2.5) and bacterial concentration was best described by a cubic polynomial. 
# While the cubic term was marginally significant (p = 0.105), it was necessary for model convergence, 
# as simpler quadratic models failed to converge properly.

simulationOutput <- simulateResiduals(fittedModel = Smoke_BacteriaModel, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)

pm25_emmeans <- emmeans(Smoke_BacteriaModel, 
                        ~MedianPM2.5_ug.m3,  
                        at = list(MedianPM2.5_ug.m3 = seq(min(smoke_bacteria$MedianPM2.5_ug.m3, na.rm = TRUE),
                                                 max(smoke_bacteria$MedianPM2.5_ug.m3, na.rm = TRUE),
                                                 length.out = 50),
                                  TotalCells_FBLBcorr = 0),
                        type = "response")

pm25_plot_data <- as.data.frame(pm25_emmeans)

pm25_values <- pm25_plot_data$MedianPM2.5_ug.m3
predictions <- pm25_plot_data$response

cubic_fit <- lm(predictions ~ poly(log(pm25_values), 3, raw = TRUE))
coeffs <- coef(cubic_fit)

# Your calculator formula
cat("TotalSpores =", coeffs[1], "+", coeffs[2], "* log(PM2.5) +", 
    coeffs[3], "* [log(PM2.5)]^2 +", coeffs[4], "* [log(PM2.5)]^3")


ggplot() +
  geom_ribbon(data = pm25_plot_data, 
              aes(x = MedianPM2.5_ug.m3, ymin = (asymp.LCL*7225)/0.03, ymax = (asymp.UCL*7225)/0.03), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = pm25_plot_data, 
            aes(x = MedianPM2.5_ug.m3, y = (response*7225)/0.03), 
            size = 1.2, color = "blue") +
  geom_point(data = smoke_bacteria, 
             aes(x = MedianPM2.5_ug.m3, y = TotalCells.m3), 
             alpha = 0.6, size = 2, color = "black") +
  labs(x = expression(PM[2.5]~(μg/m^3)),
       y = expression(Bacteria~m^-3),
       title = "Bacteria Concentration vs PM2.5 (Cubic Relationship)") +
  theme_minimal()





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









