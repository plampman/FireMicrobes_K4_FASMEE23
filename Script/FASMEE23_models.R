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
library(lme4)
library(nlme)
library(glmmTMB)
library(DHARMa)
library(mgcv)
# library(vegan)
library(ggstatsplot)


### FASMEE Spores

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam10x <- 2 ##mm
FOV10x_area <- ((FOV_diam10x/2)^2)*pi ##mm^2
FOV10x.filter <- Filter_area/FOV10x_area

FOV10x_10_area = FOV10x_area * 10

scaling = Filter_area/FOV10x_10_area

pi = pi

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam100x <- 0.2 ##mm
FOV100x_area <- ((FOV_diam100x/2)^2)*pi ##mm^2
FOV100x.filter <- Filter_area/FOV100x_area ##FoVs/filter

FOV100x_40_area = FOV100x_area * 40

scaling_1000x = Filter_area/FOV100x_40_area





SampleInfo <- read.csv('./Input_Data/FASMEE23/FASMEE_Sample_Info20250329.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(DATE_UTC = mdy(DATE_UTC, tz = 'MST7MDT'),
         DATE_UTC = ymd(DATE_UTC),
         SampleStart_UTC = ymd_hms(paste(DATE_UTC, SampleStart_UTC), tz = 'UTC'),
         SampleEnd_UTC = ymd_hms(paste(DATE_UTC, SampleEnd_UTC), tz = 'UTC'),
         DateTime_MDT = mdy_hm(DateTime_MDT, tz = 'MST7MDT'),
         DateTime_UTC =  mdy_hm(DateTime_UTC, tz = 'UTC'))

CellSamples <- SampleInfo %>%
  filter(FilterType == 'PTFE', SampleType == 'Smoke' | SampleType == 'Ambient') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+"))%>%
  filter(SampleRep == "B")

cells <- read.csv('./Input_Data/FASMEE23/Cell_Counts_FASMEE23_20250414.csv', header = T)

spores <- cells %>%
  filter(StainType == "CW/KOH") %>%
  select(-LIVECounts, -DEADCounts, -TotalCells) %>%
  rename(SampleID = 'SlideID') %>%
  mutate(SampleID = if_else(SampleType != "LabBlank" & SampleType != "FieldBlank", gsub("_", "", SampleID), SampleID),
         SampleID = if_else(SampleType == "FieldBlank", gsub("_A", "", SampleID), SampleID))

spores <- spores %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240503 & StainDate < 20240807 ~ 'A',
      StainDate >= 20240807 & StainDate < 20241119 ~ 'B',
      StainDate >= 20241119 & StainDate < 20241210 ~ 'C',
      ## There were two different TBE solutions with lab blanks from the same stain date, the first set below is TBE1 and second set is TBE2
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank1A_B_20241210", "R5B", "R4B", "R3B", "R2B", "R1B", "B3B", "B2B") ~ 'D',
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank2A_B_20241210", "B6B", "R6B", "R3B", "R7B", "R9B", "B9B", "R10B") ~ 'E', 
      StainDate >= 20250205 ~ 'F',
    )))

spores <- left_join(spores, SampleInfo, by = "SampleID")

spores <- spores %>%
  mutate(
    Volume_L = if_else(is.na(Volume_L), 0, Volume_L),
    RepVolume_L = Volume_L/2, # A and B slide replicates (S9PI vs CW/KOH)
    RepVolume_m3 = RepVolume_L/1000) %>%
  select(-SampleType.y, -Sampler.y) %>%
  rename(SampleType = 'SampleType.x', Sampler = 'Sampler.x')

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
  mutate(
    Platform = if_else(is.na(Platform), "Blank", Platform),
    SampleType = factor(SampleType),
    Platform = factor(Platform),
    SampleID = factor(SampleID),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)) %>% 
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient") | Platform == "Blank")

spores_blue <- spores_stat_test %>%
  filter(Platform == "Blue") %>%
  mutate(Sample_num = if_else(Platform == "Blue", str_extract(Sample, "\\d+"), NA_character_),
         Sample_num = as.numeric(Sample_num))

spores_blue_pa <- left_join(spores_blue, PA_stats_FASMEE23, by = c('Sample_num' = 'Sample'))

spores_blue_pa <- spores_blue_pa %>%
  mutate(
    presence = TotalSpores_LBcorr > 0,
    MedianMR_centered = scale(MedianMR, scale = F)
    )

spores_blue_pa_C <- left_join(spores_blue_pa, slim_fasmmee_C, by = c('Sample_num' = 'Sample'))

#write.csv(spores_blue_pa_C, 'FASMEE23_Spores_PA_C_20250501.csv', row.names = F)

sample_spores_blue <- spores_blue  %>%
  group_by(SampleID, Sample, Project, SampleType, RepVolume_m3, StainDate, DateCounted, StainType) %>%
  summarise(
    median_spores.FOV = median(TotalCells_LBcorr),
    mean_spores.FOV = mean(TotalCells),
    sd_spores.FOV = sd(TotalCells),
    total_median_spores.filter = median_spores.FOV*FOV100x.filter,
    total_mean_spores.filter = mean_spores.FOV*FOV100x.filter,
    Live.Dead = mean(LiveCells.TotalCells)
  ) %>% ungroup %>%
  mutate(
    total_median_spores.m3 = total_median_spores.filter/RepVolume_m3
    )
    
# Total spore model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------
ggplot(spores_blue_pa, aes(x = SmokeLevel, y = MedianMR)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw()

spores_blue_pa <- spores_blue_pa %>%
  filter(SampleID != "B1B")

unique(spores_blue_pa$SampleID)

TotalSpores_m <- glmmTMB(TotalSpores_LBcorr ~ SmokeLevel*MedianMR + offset(log_volume_offset_m3) + 
                          (1|SampleID),
                        family=tweedie(link="log"), data = spores_blue_pa, ziformula = ~0)
summary(TotalSpores_m)


TotalSpores_m <- glmmTMB(TotalSpores ~ SmokeLevel*poly(MedianMR, 2) + offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) + 
                           (1|LB_Batch:SampleID),
                         family=nbinom2(link="log"), data = spores_blue_pa, ziformula = ~0)
summary(TotalSpores_m)

simulationOutput <- simulateResiduals(fittedModel = TotalSpores_m, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testQuantiles(simulationOutput)


# Create a grid of MedianMR values covering your observed range (assuming 3.0-3.8)
mr_seq <- seq(3.0, 3.8, by = 0.2)

# Get predicted values at each combination of SmokeLevel and MedianMR
em_interact <- emmeans(TotalSpores_m, 
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

#Spore hurdle model
#-----------------------------------------------------------------------------------

presence_model <- glmmTMB(presence ~ SmokeLevel+MedianMR + 
                            (1|LB_Batch:SampleID),
                          family = binomial(link = "logit"), 
                          data = spores_blue_pa)
summary(presence_model)

positive_only <- subset(spores_blue_pa, TotalSpores_LBcorr > 0)

# Fit a model for positive values only
positive_model <- glmmTMB(TotalSpores_LBcorr ~ SmokeLevel+MedianMR + 
                            offset(log_volume_offset_m3) + 
                            (1|Sample),
                          family = Gamma(link = "log"), 
                          data = positive_only)
summary(positive_model)

simulationOutput <- simulateResiduals(fittedModel = positive_model, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
testDispersion(TotalCells_m)


# Spore smoke model
#-------------------------------------------------------------------------------------

smoke_spores_pa_C <- spores_blue_pa_C %>% 
  filter(SampleType == "Smoke") 

smoke_spores_pa_C$logPM25 <- log(smoke_spores_pa_C$MedianPM2.5_ug.m3)
smoke_spores_pa_C$logPM10 <- log(smoke_spores_pa_C$MedianPM10_ug.m3)
smoke_spores_pa_C$logPM1<- log(smoke_spores_pa_C$MedianPM1.0_ug.m3)

write.csv(smoke_spores_pa_C, './Output/K4_SmokeSpores.csv', row.names = F)


# model_logPM <- glmmTMB(TotalSpores_LBcorr ~ logPM25*MedianMCE + MedianMR +
#                          offset(log_volume_offset_m3) + (1|Sample),
#                        family=tweedie(link="log"), data = smoke_spores_pa_C)
# summary(model_logPM)

model_logPM <- glmmTMB(TotalSpores ~ logPM25 + MedianMCE + MedianMR + 
                         logPM25:MedianMCE + logPM25:MedianMR +
                         offset(log_volume_offset_m3) + offset(log1TotalSpores_LB) + (1|Sample),
                       family=nbinom2(link="log"), data = smoke_spores_pa_C)

summary(model_logPM)

model_logPM <- glmmTMB(TotalSpores_LBcorr ~ logPM25 + MedianMCE + MedianMR +
                         logPM25:MedianMCE + logPM25:MedianMR +
                         offset(log_volume_offset_m3) + (1|Sample),
                       family=tweedie(link="log"), data = smoke_spores_pa_C)

summary(model_logPM)

simulationOutput <- simulateResiduals(fittedModel = model_logPM, plot = F)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)



########### FASMEE Bacteria ###########################################################
#######################################################################################

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam100x <- 0.2 ##mm
FOV100x_area <- ((FOV_diam100x/2)^2)*pi ##mm^2
FOV100x.filter <- Filter_area/FOV100x_area ##FoVs/filter

SampleInfo <- read.csv('./Input_Data/FASMEE23/FASMEE_Sample_Info20250329.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(DATE_UTC = mdy(DATE_UTC, tz = 'MST7MDT'),
         DATE_UTC = ymd(DATE_UTC),
         SampleStart_UTC = ymd_hms(paste(DATE_UTC, SampleStart_UTC), tz = 'UTC'),
         SampleEnd_UTC = ymd_hms(paste(DATE_UTC, SampleEnd_UTC), tz = 'UTC'),
         DateTime_MDT = mdy_hm(DateTime_MDT, tz = 'MST7MDT'),
         DateTime_UTC =  mdy_hm(DateTime_UTC, tz = 'UTC'))

CellSamples <- SampleInfo %>%
  filter(FilterType == 'PTFE', SampleType == 'Smoke' | SampleType == 'Ambient') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+"))%>%
  filter(SampleRep == "B")

blue_info <- SampleInfo %>%
  filter(Platform == 'Blue',
         FilterType == 'PTFE') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+")) %>%
  filter(SampleRep == "A") %>%
  mutate(int = seq(1, by = 2, length.out = n()))

cells <- read.csv('./Input_Data/FASMEE23/Cell_Counts_FASMEE23_20250414.csv', header = T)

bacteria <- cells %>%
  filter(StainType == 'S9PI') %>%
  rename(SampleID = 'SlideID') %>%
  mutate(SampleID = if_else(SampleType != "LabBlank" & SampleType != "FieldBlank", gsub("_", "", SampleID), SampleID),
         SampleID = if_else(SampleType == "FieldBlank", gsub("_A", "", SampleID), SampleID))

bacteria <- bacteria %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240503 & StainDate < 20240807 ~ 'A',
      StainDate >= 20240807 & StainDate < 20241119 ~ 'B',
      StainDate >= 20241119 & StainDate < 20241210 ~ 'C',
      ## There were two different TBE solutions with lab blanks from the same stain date, the first set below is TBE1 and second set is TBE2
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank1A_A_20241210", "R5B", "R4B", "R3B", "R2B", "R1B", "B3B", "B2B") ~ 'D',
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank3A_A_20241210", "B6B", "R6B", "R3B", "R7B", "R9B", "B9B", "R10B") ~ 'E', 
      StainDate >= 20250205 ~ 'F',
    )))

bacteria <- left_join(bacteria, SampleInfo, by = "SampleID")

bacteria <- bacteria %>%
  mutate(
    Volume_L = if_else(is.na(Volume_L), 0, Volume_L),
    RepVolume_L = Volume_L/2, # A and B slide replicates (S9PI vs CW/KOH)
    RepVolume_m3 = RepVolume_L/1000) %>%
  select(-SampleType.y, -Sampler.y) %>%
  rename(SampleType = 'SampleType.x', Sampler = 'Sampler.x')

blank_means <- bacteria %>%
  filter(SampleType == "LabBlank") %>%
  group_by(LB_Batch) %>%
  summarise(
    TotalCells_LB = mean(TotalCells),
    DeadCells_LB = mean(DEADCounts)
    ) %>%
  ungroup

bacteria <- left_join(bacteria, blank_means, by = "LB_Batch")

bacteria <- bacteria %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    log1TotalCells_LB = log1p(TotalCells_LB),
    log1TotalCells = log1p(TotalCells),
    TotalCells_LBcorr = pmax(0, TotalCells - TotalCells_LB),
    log1TotalCells_LBcorr = log1p(TotalCells_LBcorr),
    TotalCells.filter_LBcorr = TotalCells_LBcorr*FOV100x.filter,
    DeadCells_LBcorr = pmax(0, DEADCounts - DeadCells_LB),
    LiveCells_LBcorr = pmax(0, TotalCells_LBcorr - DeadCells_LBcorr),
    DeadCells.filter_LBcorr = DeadCells_LBcorr*FOV100x.filter,
    LiveCells.filter_LBcorr = LiveCells_LBcorr*FOV100x.filter,
    TotalCells_LBcorr_m3 = TotalCells_LBcorr/RepVolume_m3,
    DeadCells_LBcorr_m3 = DeadCells_LBcorr/RepVolume_m3,
    LiveCells_LBcorr_m3 = LiveCells_LBcorr/RepVolume_m3,
    LiveCells.TotalCells = if_else(LiveCells_LBcorr_m3 > 0, LiveCells_LBcorr_m3/TotalCells_LBcorr_m3, 0),
    LiveCells.TotalCells_adj = case_when(
      LiveCells.TotalCells == 0 ~ 0.001,  
      LiveCells.TotalCells == 1 ~ 0.999,  
      TRUE ~ LiveCells.TotalCells))

bacteria_stat_test <- bacteria %>%
  mutate(
    Platform = if_else(is.na(Platform), "Blank", Platform),
    SampleType = factor(SampleType),
    Platform = factor(Platform),
    Day = factor(Day),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)) %>% 
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient") | Platform == "Blank")

bacteria_blue <- bacteria_stat_test %>%
  filter(Platform == "Blue") %>%
  mutate(Sample_num = if_else(Platform == "Blue", str_extract(Sample, "\\d+"), NA_character_),
         Sample_num = as.numeric(Sample_num))

bacteria_blue_pa <- left_join(bacteria_blue, PA_stats_FASMEE23, by = c('Sample_num' = 'Sample'))

bacteria_blue_pa <- bacteria_blue_pa %>%
  mutate(
    presence = TotalCells_LBcorr > 0)

#write.csv(bacteria_blue_pa, 'bacteria_blue_pa_FASMEE20250418.csv', row.names = F)
    

sample_bacteria_blue <- bacteria_blue  %>%
  group_by(SampleID, Sample, Project, SampleType, RepVolume_m3, StainDate, DateCounted, StainType) %>%
  summarise(
    median_bacteria.FOV = median(TotalCells_LBcorr),
    mean_bacteria.FOV = mean(TotalCells),
    sd_bacteria.FOV = sd(TotalCells),
    total_median_bacteria.filter = median_bacteria.FOV*FOV100x.filter,
    total_mean_bacteria.filter = mean_bacteria.FOV*FOV100x.filter,
    Live.Dead = mean(LiveCells.TotalCells)
  ) %>% ungroup %>%
  mutate(
    total_median_bacteria.m3 = total_median_bacteria.filter/RepVolume_m3
  )

# Total cell model for ambient versus smoke with mixing ratio
#----------------------------------------------------------------------------------------------------

# TotalCells_m <- glmmTMB(TotalCells_LBcorr ~ SampleType*MedianMR + offset(log_volume_offset_m3) + 
#                         (1|LB_Batch:SampleID),
#                       family=tweedie(link="log"), data = bacteria_blue_pa, ziformula = ~SampleType)
# summary(TotalCells_m)


TotalCells_m <- glmmTMB(TotalCells ~ SmokeLevel+MedianMR + 
                          offset(log_volume_offset_m3) + offset(log1TotalCells_LB) +
                          (1|SampleID),
                        family=nbinom2(link="log"), dispformula = ~SmokeLevel, 
                        data = bacteria_blue_pa, ziformula = ~0)
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

