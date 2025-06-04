#########################################################################
###
### Program name: Konza 4 spore analysis
###
### Purpose: Konza 4 spore concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 01/29/2025
###
#########################################################################

library(tidyverse)

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam1000x <- 0.2 ##mm
FOV1000x_area <- ((FOV_diam1000x/2)^2)*pi ##mm^2
FOV1000x.filter <- Filter_area/FOV1000x_area ##FoVs/filter

FOV1000x_40_area = FOV1000x_area * 40
scaling_1000x = Filter_area/FOV1000x_40_area


spores <- read_csv('./Input_Data/Konza4/Konza4_SporeCounts_1000X_20250513.csv') %>%
  left_join(., SlideInfo, by = "SlideID") %>%
  mutate(SampleID = if_else(SampleType != "LabBlank", gsub("_.*$", "", SlideID), SlideID))

spores <- spores %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240621 & StainDate < 20240703 ~ 'K4_A',
      StainDate >= 20240703 & StainDate < 20241119 ~ 'K4_B',
      StainDate >= 20241119 & StainDate < 20250329 ~ 'K4_C',
      StainDate >= 20250329 ~ 'K4_D',
      TRUE ~ NA_character_
    ))
  )

spores <- left_join(spores, SampleInfo, by = c("SampleID" = "Sample"))

spores <- spores %>%
  mutate(
    RepVolume_m3 = Slide_RepVolume_L/1000)

blank_means <- spores %>%
  filter(SampleType == "LabBlank") %>%
  group_by(LB_Batch) %>%
  summarise(
    TotalSpores_LB = mean(TotalSpores),
    TotalSpores_LB.filter = TotalSpores_LB*FOV1000x.filter,
    log1TotalSpores_LB = log1p(TotalSpores_LB),
    log1TotalSpores_LB.filter = log1p(TotalSpores_LB.filter) 
  ) %>%
  ungroup

spores <- left_join(spores, blank_means, by = "LB_Batch")

spores <- spores %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    TotalSpores_LBcorr = pmax(0, TotalSpores - TotalSpores_LB),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)) 

FieldBlank_mean <- spores %>%
  filter(SampleType == "FieldBlank") %>%
  summarise(
    TotalSpores_FB = mean(TotalSpores_LBcorr)
  ) %>%
  ungroup

spores <- spores %>%
  filter(SampleType != "FieldBlank") %>%
  mutate(
    TotalSpores_FBLBcorr = pmax(0, TotalSpores_LBcorr - FieldBlank_mean$TotalSpores_FB))

Ambient_mean <- spores %>%
  filter(SampleType == "Ambient") %>%
  summarise( # Ambient mean for background correction
    AmbientSpores.FOV = mean(TotalSpores_FBLBcorr)
  )

spores <- spores %>%
  filter(SampleType != "FieldBlank") %>%
  mutate(Sample_num = as.numeric(str_extract(SampleID, "\\d+")),
         TotalSpores.filter = TotalSpores*FOV1000x.filter,
         TotalSpores.filter_FBLBcorr = TotalSpores_FBLBcorr*FOV1000x.filter,
         TotalSpores_FBLBcorr.m3 = TotalSpores.filter_FBLBcorr/RepVolume_m3,
         TotalSpores_Bcorr = if_else(SampleType == "Smoke", pmax(0, TotalSpores_FBLBcorr - Ambient_mean$AmbientSpores.FOV), NA),
         TotalSpores_Bcorr.m3 = (TotalSpores_Bcorr*FOV1000x.filter)/RepVolume_m3,
         log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)
  )

spores_pa <- left_join(spores, PA_stats_k4, by = c('Sample_num' = 'Sample'))

spores_pa_C <- left_join(spores_pa, slim_UI_EPA_C, by = c('Sample_num' = 'Sample')) %>%
  mutate(
    spores.Mg = TotalSpores_Bcorr.m3/biomass_Mg,
    log1spores.Mg = log1p(spores.Mg)
  )

na_count <- spores_pa_C %>%
  summarize(across(everything(), ~sum(is.na(.))))

#write.csv(spores_pa_C, './Output/Output_data/K4/k4_Spores_PA_C.csv', row.names = F)

sample_spores <- spores_pa_C  %>%
  group_by(SampleID, Unit, AQI_PM2.5, RepVolume_m3) %>%
  summarise(
    mean_Spores.Mg = mean(spores.Mg, na.rm = T),
    meanMCE = mean(MeanMCE, na.rm = T),
    meanMR = mean(MedianMR),
    mean_spores.FOV = mean(TotalSpores),
    sd_spores.FOV = sd(TotalSpores),
    mean_spores.m3 = mean(TotalSpores_FBLBcorr.m3),
    mean_BcorrSpores.m3 = mean(TotalSpores_Bcorr.m3, na.rm = T),
    sd_spores.m3 = sd(TotalSpores_FBLBcorr.m3))

FEF_K4 <- spores_pa_C  %>%
  summarise(
    Spores.Mg = mean(spores.Mg, na.rm = T),
    Spores.m3 = mean(TotalSpores_FBLBcorr.m3))






