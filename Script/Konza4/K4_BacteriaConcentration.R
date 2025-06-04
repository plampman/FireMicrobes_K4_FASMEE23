#########################################################################
###
### Program name: Konza 4 bacteria analysis
###
### Purpose: Konza 4 bacteria concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/17/2025
###
#########################################################################

library(tidyverse)
library(gridExtra)
library(ggstatsplot)

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam1000x <- 0.2 ##mm
FOV1000x_area <- ((FOV_diam1000x/2)^2)*pi ##mm^2
FOV1000x.filter <- Filter_area/FOV1000x_area ##FoVs/filter

FOV1000x_40_area = FOV1000x_area * 40
scaling_1000x = Filter_area/FOV1000x_40_area

bacteria <- read.csv('./Input_Data/Konza4/K4_BacteriaCounts_1000X_20250401.csv', header = T)
amb_cells_k3 <- read.csv('./Input_Data/Konza4/Konza3_Bacteria_Ambient_FBLBcorr.csv', header = T)

bacteria <- bacteria %>%
  filter(SlideID != "10B_A") %>% #Issues with this sample due to staining--REMOVE ALWAYS
  mutate(SampleID = if_else(SampleType != "LabBlank", gsub("_.*$", "", SlideID), SlideID))

bacteria <- bacteria %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240619 & StainDate < 20240703 ~ 'K4_A',
      StainDate >= 20240703 & StainDate < 20241119 ~ 'K4_B',
      StainDate >= 20241119 & StainDate < 20250329 ~ 'K4_C', # This only applies to K4 spores
      StainDate >= 20250329 ~ 'K4_D',
      TRUE ~ NA_character_
    )))

bacteria <- left_join(bacteria, SampleInfo, by = c("SampleID" = "Sample"))

bacteria <- bacteria %>%
  mutate(
    RepVolume_m3 = Slide_RepVolume_L/1000
    )

blank_means <- bacteria %>%
  filter(SampleType == "LabBlank") %>%
  group_by(LB_Batch) %>%
  summarise(
    TotalCells_LB = mean(TotalCells),
    LiveCells_LB = mean(LIVECounts),
    DeadCells_LB = mean(DEADCounts)
  ) %>%
  ungroup

bacteria <- left_join(bacteria, blank_means, by = "LB_Batch")

bacteria <- bacteria %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    DeadCells_LBcorr = pmax(0, DEADCounts - DeadCells_LB),
    LiveCells_LBcorr = pmax(0, LIVECounts - LiveCells_LB),
    TotalCells_LBcorr = pmax(0, TotalCells - TotalCells_LB))

FieldBlank_mean <- bacteria %>%
  filter(SampleType == "FieldBlank") %>%
  summarise(
    TotalCells_FB = mean(TotalCells_LBcorr),
    LiveCells_FB = mean(LiveCells_LBcorr),
    DeadCells_FB = mean(DeadCells_LBcorr)
  ) %>%
  ungroup

bacteria <- bacteria %>%
  filter(SampleType != "FieldBlank") %>%
  mutate(
    DeadCells_FBLBcorr = pmax(0, DeadCells_LBcorr - FieldBlank_mean$DeadCells_FB),
    LiveCells_FBLBcorr = pmax(0, LiveCells_LBcorr - FieldBlank_mean$LiveCells_FB),
    TotalCells_FBLBcorr = pmax(0, TotalCells_LBcorr - FieldBlank_mean$TotalCells_FB),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0))


Ambient_mean <- amb_cells_k3 %>% # Using Konza 3 ambient samples for background corrections
  summarise( # Ambient mean for background correction
    AmbientTotalCells.FOV = mean(Totalbacteria.FOV_FBLBcorr),
    AmbientLiveCells.FOV = mean(Livebacteria.FOV_FBLBcorr),
    AmbientDeadCells.FOV = mean(Deadbacteria.FOV_FBLBcorr)
  )

bacteria <- bacteria %>%
  mutate(
    TotalCells_Bcorr = if_else(SampleType == "Smoke", pmax(0, TotalCells_FBLBcorr - Ambient_mean$AmbientTotalCells.FOV), NA),
    LiveCells_Bcorr = if_else(SampleType == "Smoke", pmax(0, LiveCells_FBLBcorr - Ambient_mean$AmbientLiveCells.FOV), NA),
    DeadCells_Bcorr = if_else(SampleType == "Smoke", pmax(0, DeadCells_FBLBcorr - Ambient_mean$AmbientDeadCells.FOV), NA),
    TotalCells.filter_Bcorr = TotalCells_Bcorr*FOV1000x.filter,
    DeadCells.filter_Bcorr = DeadCells_Bcorr*FOV1000x.filter,
    LiveCells.filter_Bcorr = LiveCells_Bcorr*FOV1000x.filter,
    TotalCells_Bcorr.m3 = TotalCells.filter_Bcorr/RepVolume_m3,
    DeadCells_Bcorr.m3 = DeadCells.filter_Bcorr/RepVolume_m3,
    LiveCells_Bcorr.m3 = LiveCells.filter_Bcorr/RepVolume_m3,
    
    Live.Total_Bcorr = if_else(LiveCells_Bcorr.m3 > 0, LiveCells_Bcorr.m3/TotalCells_Bcorr.m3, 0),
    
    TotalCells.filter_FBLBcorr = TotalCells_FBLBcorr*FOV1000x.filter,
    DeadCells.filter_FBLBcorr = DeadCells_FBLBcorr*FOV1000x.filter,
    LiveCells.filter_FBLBcorr = LiveCells_FBLBcorr*FOV1000x.filter,
    TotalCells_FBLBcorr.m3 = TotalCells.filter_FBLBcorr/RepVolume_m3,
    DeadCells_FBLBcorr.m3 = DeadCells.filter_FBLBcorr/RepVolume_m3,
    LiveCells_FBLBcorr.m3 = LiveCells.filter_FBLBcorr/RepVolume_m3,
    
    LiveCells.TotalCells = if_else(LiveCells_FBLBcorr.m3 > 0, LiveCells_FBLBcorr.m3/TotalCells_FBLBcorr.m3, 0)
  )

bacteria <- bacteria %>%
  mutate(
    SampleID = factor(SampleID),
    SampleType = factor(SampleType),
    Day = factor(Day),
    log_volume_offset_m3 = log(RepVolume_m3))

bacteria_pa <- left_join(bacteria, PA_stats_k4, by = c('Sample_num' = 'Sample'))

bacteria_pa_C <- left_join(bacteria_pa, slim_UI_EPA_C, by = c('Sample_num' = 'Sample')) %>%
  mutate(
    Total_bacteria.Mg = TotalCells_Bcorr.m3/biomass_Mg,
    Live_bacteria.Mg = LiveCells_Bcorr.m3/biomass_Mg
  )

write.csv(bacteria_pa_C, './Output/Output_data/K4/k4_Bacteria_PA_C.csv', row.names = F)

na_count <- bacteria_pa_C %>%
  summarize(across(everything(), ~sum(is.na(.))))

sample_bacteria <- bacteria_pa_C  %>%
  group_by(SampleID, Unit, AQI_PM2.5, RepVolume_m3) %>%
  summarise(
    
    meanMCE = mean(MeanMCE, na.rm = T),
    meanMR = mean(MedianMR),
    
    TotalBacteria.Mg = mean(Total_bacteria.Mg, na.rm = T),
    LiveBacteria.Mg = mean(Live_bacteria.Mg, na.rm = T),
    Total_BcorrBacteria.m3 = mean(TotalCells_Bcorr.m3),
    Live_BcorrBacteria.m3 = mean(LiveCells_Bcorr.m3),
    sd_total_BcorrBacteria.m3 = sd(TotalCells_Bcorr.m3),
    sd_live_BcorrBacteria.m3 = sd(LiveCells_Bcorr.m3),
    Total_Bacteria.m3 = mean(TotalCells_FBLBcorr.m3),
    Live_bacteria.m3 = mean(LiveCells_FBLBcorr.m3),
    sd_total_bacteria.m3 = sd(TotalCells_FBLBcorr.m3),
    sd_live_bacteria.m3 = sd(LiveCells_FBLBcorr.m3))


BEF_K4 <- bacteria_pa_C  %>%
  summarise(
    TotalBacteria.Mg = mean(Total_bacteria.Mg, na.rm = T),
    TotalBacteria.m3 = mean(TotalCells_FBLBcorr.m3))

    
    
    