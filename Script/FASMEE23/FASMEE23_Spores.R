#########################################################################
###
### Program name: FASMEE23 spore analysis
###
### Purpose: FASMEE23 spore concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/25/2025
###
#########################################################################

# Decimals in names "." denote "/" 

#Spore counts from Liv Lampman

library(tidyverse)


Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam1000x <- 0.2 ##mm
FOV1000x_area <- ((FOV_diam1000x/2)^2)*pi ##mm^2
FOV1000x.filter <- Filter_area/FOV1000x_area ##FoVs/filter

FOV1000x_40_area = FOV1000x_area * 40
scaling_1000x = Filter_area/FOV1000x_40_area


spores <- read_csv('./Input_Data/FASMEE23/FASMEE23_SporeCounts_1000X_20250513.csv') %>%
  left_join(., SlideInfo, by = "SlideID") %>%
  mutate(SampleID = if_else(SampleType != "LabBlank", gsub("_.*$", "", SlideID), SlideID))

spores <- spores %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240503 & StainDate < 20240807 ~ 'F23_A',
      StainDate >= 20240807 & StainDate < 20241119 ~ 'F23_B',
      StainDate >= 20241119 & StainDate < 20241210 ~ 'F23_C',
      ## There were two different TBE solutions with lab blanks from the same stain date, the first set below is TBE1 and second set is TBE2
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank1A_B_20241210", "R5B", "R4B", "R3B", "R2B", "R1B", "B3B", "B2B") ~ 'F23_D',
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank2A_B_20241210", "B6B", "R6B", "R3B", "R7B", "R9B", "B9B", "R10B") ~ 'F23_E', 
      StainDate >= 20250205 ~ 'F23_F',
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
    TotalSpores_LB = mean(TotalSpores),
    log1TotalSpores_LB = log1p(TotalSpores_LB),
  ) %>%
  ungroup


spores <- left_join(spores, blank_means, by = "LB_Batch")

spores <- spores %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    TotalSpores_LBcorr = pmax(0, TotalSpores - TotalSpores_LB),
    TotalSpores_FBLBcorr = pmax(0, TotalSpores_LBcorr - FieldBlank_mean$TotalSpores_FB),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)) %>% 
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient") | Platform == "Blank")

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
  mutate(Sample_num = if_else(Platform == "Blue", str_extract(SampleID, "\\d+"), NA_character_),
         Sample_num = as.numeric(Sample_num),
         TotalSpores.filter = TotalSpores*FOV1000x.filter,
         TotalSpores.filter_FBLBcorr = TotalSpores_FBLBcorr*FOV1000x.filter,
         TotalSpores_FBLBcorr.m3 = TotalSpores.filter_FBLBcorr/RepVolume_m3,
         TotalSpores_Bcorr = if_else(SampleType == "Smoke", pmax(0, TotalSpores_FBLBcorr - Ambient_mean$AmbientSpores.FOV), NA),
         TotalSpores_Bcorr.m3 = (TotalSpores_Bcorr*FOV1000x.filter)/RepVolume_m3,
         log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)
         )

spores_pa <- left_join(spores, PA_stats_FASMEE23, by = c('Sample_num' = 'Sample'))

spores_pa_C <- left_join(spores_pa, slim_fasmmee_C, by = c('Sample_num' = 'Sample')) %>%
  mutate(
    Platform = if_else(is.na(Platform), "Blank", Platform),
    SampleType = factor(SampleType),
    Platform = factor(Platform),
    Day = factor(Day),
    SampleID = factor(SampleID),
    spores.kg = TotalSpores_Bcorr.m3/biomass_kg)

unique(spores_pa_C$SampleID)

na_count <- spores_pa_C %>%
  summarize(across(everything(), ~sum(is.na(.))))

#write.csv(spores_blue_pa_C, 'FASMEE23_Spores_PA_C_20250501.csv', row.names = F)

sample_spores <- spores_pa_C  %>%
  group_by(SampleID, AQI_PM2.5, RepVolume_m3) %>%
  summarise(
    meanMCE = mean(MeanMCE, na.rm = T),
    meanlogPM2.5 = mean(logPM2.5),
    meanMR = mean(MedianMR),
    mean_spores.FOV = mean(TotalSpores),
    mean_sporesFBLBcorr.FOV = mean(TotalSpores_FBLBcorr),
    sd_spores.FOV = sd(TotalSpores),
    mean_spores.m3 = mean(TotalSpores_FBLBcorr.m3),
    mean_BcorrSpores.m3 = mean(TotalSpores_Bcorr.m3, na.rm = T),
    mean_Spores.kg = mean(spores.kg, na.rm = T),
    sd_spores.m3 = sd(TotalSpores_FBLBcorr.m3))

