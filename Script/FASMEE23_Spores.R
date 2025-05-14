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
library(gridExtra)
library(ggstatsplot)


Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam100x <- 0.2 ##mm
FOV100x_area <- ((FOV_diam100x/2)^2)*pi ##mm^2
FOV100x.filter <- Filter_area/FOV100x_area ##FoVs/filter

FOV100x_40_area = FOV100x_area * 40
scaling_1000x = Filter_area/FOV100x_40_area


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
    TotalSpores_LB = mean(TotalSpores)
  ) %>%
  ungroup

spores <- left_join(spores, blank_means, by = "LB_Batch")

spores <- spores %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    log1TotalSpores_LB = log1p(TotalSpores_LB),
    TotalSpores_LBcorr = pmax(0, TotalSpores - TotalSpores_LB),
    TotalSpores.filter_LBcorr = TotalSpores_LBcorr*FOV100x.filter,
    TotalSpores_LBcorr_m3 = TotalSpores_LBcorr/RepVolume_m3,
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)) %>% 
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient") | Platform == "Blank")

spores_blue <- spores %>%
  filter(Platform == "Blue") %>%
  filter(SampleType != "FieldBlank") %>%
  mutate(Sample_num = if_else(Platform == "Blue", str_extract(SampleID, "\\d+"), NA_character_),
         Sample_num = as.numeric(Sample_num))

spores_blue_pa <- left_join(spores_blue, PA_stats_FASMEE23, by = c('Sample_num' = 'Sample'))

spores_blue_pa_C <- left_join(spores_blue_pa, slim_fasmmee_C, by = c('Sample_num' = 'Sample')) %>%
  mutate(
    Platform = if_else(is.na(Platform), "Blank", Platform),
    SampleType = factor(SampleType),
    Platform = factor(Platform),
    SampleID = factor(SampleID))



unique(spores_blue_pa_C$SampleID)

#write.csv(spores_blue_pa_C, 'FASMEE23_Spores_PA_C_20250501.csv', row.names = F)

sample_spores_blue <- spores_blue_pa_C  %>%
  group_by(SampleID, SmokeLevel, RepVolume_m3, StainDate, DateCounted, StainType) %>%
  summarise(
    median_spores.FOV = median(TotalSpores),
    mean_spores.FOV = mean(TotalSpores),
    sd_spores.FOV = sd(TotalSpores),
    total_median_spores.filter = median_spores.FOV*FOV100x.filter,
    total_mean_spores.filter = mean_spores.FOV*FOV100x.filter) %>% 
  ungroup %>%
  mutate(
    total_median_spores.m3 = total_median_spores.filter/RepVolume_m3
  )


