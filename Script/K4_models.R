#########################################################################
###
### Program name: Konza 4 statistics
###
### Purpose: Konza 4 mixed models
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 04/06/2025
###
#########################################################################

library(tidyverse)
library(emmeans)
library(lme4)
library(nlme)
library(glmmTMB)
library(DHARMa)
# library(mgcv)
# library(vegan)
library(ggstatsplot)

#Konza 4 Fungi
#--------------------------------------------------------------------------------------

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam10x <- 2 ##mm
FOV10x_area <- ((FOV_diam10x/2)^2)*pi ##mm^2
FOV10x.filter <- Filter_area/FOV10x_area

SampleInfo <- read.csv('./K4_SampleInfo.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(Sample = as.character(Sample),
         Sample_num = as.numeric(str_extract(Sample, "\\d+")))

cells <- read.csv('./K4_CELL_COUNTS_20250401.csv', header = T)

spores <- cells %>%
  filter(StainType == "CW/KOH") %>%
  select(-LIVECounts, -DEADCounts, -TotalCells, -Clumps.4ormorecells.) %>%
  mutate(Sample = as.character(Sample),
         Sample = str_trim(Sample),
         Sample_num = if_else(SampleType == "Smoke" | SampleType == "Ambient", 
                              as.numeric(str_extract(Sample, "\\d+")), NA_real_))

spores <- left_join(spores, SampleInfo, by = "Sample_num")

spores_plt <- spores %>% filter(SampleType == "Smoke" | SampleType == "Ambient")
ggplot(spores_plt, aes(x = TotalSpores, fill = SampleType)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black") +
  labs(title = "Histogram of Total Spores per FOV by Sample Type",
       x = "Total Spores/FOV",
       y = "Frequency") +
  scale_fill_manual(values = c("Ambient" = "darkgreen", "Smoke" = "orange")) +
  theme_minimal()

spores <- spores %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240621 & StainDate < 20240703 ~ 'A',
      (StainDate >= 20240703 & StainDate < 20241119 & SlideID %in% c("LabBlankE", "12A", "13A", "14A")) ~ 'B',
      (StainDate >= 20240703 & StainDate < 20241119 & !(SlideID %in% c("LabBlankE", "12A", "13A", "14A"))) ~ 'C',
      StainDate >= 20241119 & StainDate < 20250329 ~ 'D',
      StainDate >= 20250329 ~ 'E',
      TRUE ~ NA_character_
    ))
  )

spores <- spores %>%
  mutate(
    RepVolume_m3 = Slide_RepVolume_L/1000
    ) %>%
  select(-Sample.y) %>%
  rename(Sample = 'Sample.x')

blank_means <- spores %>%
  filter(SampleType == "LabBlank") %>%
  group_by(LB_Batch) %>%
  summarise(
    TotalSpores_LB = mean(TotalSpores)
  ) %>%
  ungroup

spores <- left_join(spores, blank_means, by = "LB_Batch")

spores <- spores %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    log1TotalSpores_LB = log1p(TotalSpores_LB),
    TotalSpores_LBcorr = pmax(0, TotalSpores - TotalSpores_LB),
    TotalSpores.filter_LBcorr = TotalSpores_LBcorr*FOV10x.filter,
    TotalSpores_LBcorr_m3 = TotalSpores_LBcorr/RepVolume_m3)

spores_stat_test <- spores %>%
  filter(SampleType == "Ambient" | SampleType == "Smoke") %>%
  filter(Sample != '14A') %>%
  mutate(
    SampleType = factor(SampleType),
    Unit = factor(Unit),
    Day = factor(Day),
    log_volume_offset_m3 = log(RepVolume_m3))

spores_pa <- left_join(spores_stat_test, PA_stats_k4, by = c('Sample_num' = 'Sample'))

spores_pa_C <- left_join(spores_pa, slim_UI_EPA_C, by = c('Sample_num' = 'Sample'))

# spores_pa <- spores_pa %>%
#   filter(SlideID != "14A")

spores_pa_C <- spores_pa_C %>%
  mutate(
    presence = TotalSpores_LBcorr > 0)

#write.csv(spores_pa_C, './Konza4_spores_PA_C_20250501.csv', row.names = F)


sample_spores <- spores_pa_C %>%
  group_by(Sample, Sample_num, Project, SampleType, RepVolume_m3, StainDate, DateCounted, StainType) %>%
  summarise(
    median_spores.FOV = median(TotalSpores_LBcorr),
    mean_spores.FOV = mean(TotalSpores),
    sd_spores.FOV = sd(TotalSpores),
    total_median_spores.filter = median_spores.FOV*FOV10x.filter,
    total_mean_spores.filter = mean_spores.FOV*FOV10x.filter,
  ) %>% ungroup %>%
  mutate(
    total_median_spores.m3 = total_median_spores.filter/RepVolume_m3)

# Total cell model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------


# TotalSpores_m <- glmmTMB(TotalSpores_LBcorr ~ SmokeLevel + MedianMR +
#                            offset(log_volume_offset_m3) + (1|Sample),
#                          family=tweedie(link="log"), data = spores_pa_C, ziformula = ~0)
# summary(TotalSpores_m)

TotalSpores_m <- glmmTMB(TotalSpores ~ SmokeLevel + MedianMR + 
                           offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) + (1|Sample),
                         family=nbinom2(link="log"), dispformula = ~SmokeLevel, data = spores_pa_C, ziformula = ~0)
summary(TotalSpores_m)

simulationOutput <- simulateResiduals(fittedModel = TotalSpores_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


em <- emmeans(TotalSpores_m, ~ SmokeLevel, type = "response")
# View the predicted means
summary(em, infer = T)

# Convert to ratios instead of differences
contrasts_ratio <- contrast(em, method = "trt.vs.ctrl", ref = "None", type = "response")
summary(contrasts_ratio, infer = TRUE)

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
