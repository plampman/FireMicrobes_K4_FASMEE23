#########################################################################
###
### Program name: FASMEE23 bacteria analysis
###
### Purpose: FASMEE23 bacteria concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/25/2025
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

bacteria <- read.csv('./Input_Data/FASMEE23/FASMEE23_BacteriaCounts_1000X_20250414.csv', header = T)

bacteria <- bacteria %>%
  rename(SampleID = 'SlideID') %>%
  mutate(SampleID = if_else(SampleType != "LabBlank" & SampleType != "FieldBlank", gsub("_", "", SampleID), SampleID),
         SampleID = if_else(SampleType == "FieldBlank", gsub("_A", "", SampleID), SampleID))

bacteria <- bacteria %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240503 & StainDate < 20240807 ~ 'F23A',
      StainDate >= 20240807 & StainDate < 20241119 ~ 'F23B',
      StainDate >= 20241119 & StainDate < 20241210 ~ 'F23C',
      ## There were two different TBE solutions with lab blanks from the same stain date, the first set below is TBE1 and second set is TBE2
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank1A_A_20241210", "R5B", "R4B", "R3B", "R2B", "R1B", "B3B", "B2B") ~ 'F23D',
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank3A_A_20241210", "B6B", "R6B", "R3B", "R7B", "R9B", "B9B", "R10B") ~ 'F23E', 
      StainDate >= 20250205 ~ 'F23F',
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
    DeadCells_LBcorr = pmax(0, DEADCounts - DeadCells_LB),
    TotalCells_LBcorr = pmax(0, TotalCells - TotalCells_LB))
    
FieldBlank_mean <- bacteria %>%
  filter(SampleType == "FieldBlank") %>%
  summarise(
    TotalCells_FB = mean(TotalCells_LBcorr),
    DeadCells_FB = mean(DeadCells_LBcorr)
  ) %>%
  ungroup

bacteria <- bacteria %>%
  filter(SampleType != "FieldBlank") %>%
  mutate(
    DeadCells_FBLBcorr = pmax(0, DeadCells_LBcorr - FieldBlank_mean$DeadCells_FB),
    TotalCells_FBLBcorr = pmax(0, TotalCells_LBcorr - FieldBlank_mean$TotalCells_FB),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0))


Ambient_mean <- bacteria %>%
  filter(SampleType == "Ambient") %>%
  summarise( # Ambient mean for background correction
    AmbientTotalCells.FOV = mean(TotalCells_FBLBcorr),
    AmbientDeadCells.FOV = mean(DeadCells_FBLBcorr)
  )

bacteria <- bacteria %>%
  mutate(
    TotalCells_Bcorr = if_else(SampleType == "Smoke", pmax(0, TotalCells_FBLBcorr - Ambient_mean$AmbientTotalCells.FOV), NA),
    DeadCells_Bcorr = if_else(SampleType == "Smoke", pmax(0, TotalCells_FBLBcorr - Ambient_mean$AmbientTotalCells.FOV), NA),
    TotalCells.filter_Bcorr = TotalCells_Bcorr*FOV1000x.filter,
    DeadCells.filter_Bcorr = DeadCells_Bcorr*FOV1000x.filter,
    LiveCells.filter_Bcorr = if_else(SampleType == "Smoke", TotalCells.filter_Bcorr - DeadCells.filter_Bcorr, NA),
    TotalCells_Bcorr.m3 = TotalCells.filter_Bcorr/RepVolume_m3,
    DeadCells_Bcorr.m3 = DeadCells.filter_Bcorr/RepVolume_m3,
    LiveCells_Bcorr.m3 = LiveCells.filter_Bcorr/RepVolume_m3,
    
    TotalCells.filter_FBLBcorr = TotalCells_FBLBcorr*FOV1000x.filter,
    DeadCells.filter_FBLBcorr = DeadCells_FBLBcorr*FOV1000x.filter,
    LiveCells.filter_FBLBcorr = TotalCells.filter_FBLBcorr - DeadCells.filter_FBLBcorr,
    TotalCells_FBLBcorr.m3 = TotalCells.filter_FBLBcorr/RepVolume_m3,
    DeadCells_FBLBcorr.m3 = DeadCells.filter_FBLBcorr/RepVolume_m3,
    LiveCells_FBLBcorr.m3 = LiveCells.filter_FBLBcorr/RepVolume_m3,
    
    LiveCells.TotalCells = if_else(LiveCells_FBLBcorr.m3 > 0, LiveCells_FBLBcorr.m3/TotalCells_FBLBcorr.m3, 0)
    )

bacteria_stat_test <- bacteria %>%
  mutate(
    SampleID = factor(SampleID),
    SampleType = factor(SampleType),
    Platform = factor(Platform),
    Day = factor(Day),
    log_volume_offset_m3 = log(RepVolume_m3)) %>% 
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient") | Platform == "Blank")

bacteria_blue <- bacteria_stat_test %>%
  filter(Platform == "Blue") %>%
  mutate(Sample_num = if_else(Platform == "Blue", str_extract(Sample, "\\d+"), NA_character_),
         Sample_num = as.numeric(Sample_num))

bacteria_blue_pa <- left_join(bacteria_blue, PA_stats_FASMEE23, by = c('Sample_num' = 'Sample'))

bacteria_blue_pa_C <- left_join(bacteria_blue_pa, slim_fasmmee_C, by = c('Sample_num' = 'Sample'))

na_count <- bacteria_blue_pa_C %>%
  summarize(across(everything(), ~sum(is.na(.))))

sample_bacteria_blue <- bacteria_blue_pa_C  %>%
  group_by(SampleID, SampleType, RepVolume_m3) %>%
  summarise(
    RH = mean(MeanRH, na.rm = T),
    PM2.5 = mean(MedianPM2.5_ug.m3, na.rm = T),
    Temp = mean(MeanTemp_C, na.rm = T),
    FOVs = n(),
    NAFOVs = sum(is.na(FOV)),
    meanMCE = mean(MeanMCE, na.rm = T),
    meanlogPM2.5 = mean(logPM2.5, na.rm = T),
    median_bacteria.FOV = median(TotalCells_FBLBcorr),
    mean_bacteria.FOV = mean(TotalCells_FBLBcorr),
    sd_bacteria.FOV = sd(TotalCells_FBLBcorr),
    total_mean_bacteria.filter = mean(TotalCells.filter_FBLBcorr),
    total_mean_bacteria.m3 = mean(TotalCells_FBLBcorr.m3),
    total_mean_bacteria.m3_bcorr = mean(TotalCells_Bcorr.m3),
    Live.Dead = mean(LiveCells.TotalCells)
  ) %>% ungroup




