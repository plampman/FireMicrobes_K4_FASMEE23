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
library(gridExtra)
library(ggstatsplot)

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

sample_bacteria <- bacteria %>%
  filter(Platform == "Blue") %>%
  group_by(SampleID, Sample, Project, SampleType, StainDate, DateCounted, StainType) %>%
  summarise(
    median_bacteria.FOV = median(TotalCells),
    mean_bacteria.FOV = mean(TotalCells),
    sd_bacteria.FOV = sd(TotalCells),
    total_median_bacteria.filter = median_bacteria.FOV*FOV100x.filter,
    total_mean_bacteria.filter = mean_bacteria.FOV*FOV100x.filter
  ) %>% ungroup

# sample_bacteria <- left_join(sample_bacteria, SampleInfo, by = c("SampleID", "SampleType"))
# 
# sample_bacteria <- sample_bacteria %>%
#   mutate(
#     RepVolume_L = Volume_L/2, # A and B slide replicates (S9PI vs CW/KOH)
#     RepVolume_m3 = RepVolume_L/1000,
#     Total_bacteria.L = total_bacteria.filter/RepVolume_L,
#     Total_bacteria.m3 = total_bacteria.filter/RepVolume_m3)

bacteria_LabBlanks <- sample_bacteria %>%
  filter(SampleType == "LabBlank") %>%
  group_by(StainDate) %>%
  summarise(
    "LabBlank_TotAvg" = mean(total_mean_bacteria.filter)) %>%
  mutate(StainDate = paste("Blank", StainDate, sep = "_")) %>%
  pivot_wider(names_from = StainDate, values_from = LabBlank_TotAvg)

bacteria_LBcorr <- sample_bacteria %>%
  filter(SampleType != "LabBlank") %>% 
  rowwise() %>%
  mutate(
    Totalbacteria.filter_LBcorr = as.integer(case_when(
      StainDate >= 20240503 & StainDate < 20240807 ~ total_mean_bacteria.filter - bacteria_LabBlanks$Blank_20240503,
      StainDate >= 20240807 & StainDate < 20241119 ~ total_mean_bacteria.filter - bacteria_LabBlanks$Blank_20240807,
      StainDate >= 20241119 & StainDate < 20241210 ~ total_mean_bacteria.filter - bacteria_LabBlanks$Blank_20241119,
      ## There were two different TBE solutions with lab blanks from the same stain date, the first set below is TBE1 and second set is TBE2
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% c("R5B", "R4B", "R3B", "R2B", "R1B", "B3B", "B2B") ~ total_mean_bacteria.filter - 0,
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% c("B6B", "R6B", "R3B", "R7B", "R9B", "B9B", "R10B") ~ total_mean_bacteria.filter - 7225,
      
      StainDate >= 20250205 ~ total_mean_bacteria.filter - bacteria_LabBlanks$Blank_20250205,
    ))) %>%
  ungroup

bacteria_FieldBlanks <- bacteria_LBcorr %>%
  filter(SampleType == "FieldBlank") %>%
  #group_by(DATE_UTC) %>%
  summarise(
    "FieldBlank_TotAvg" = mean(Tot_bacteria.filter))

bacteria_LB_FB_corr <- bacteria_LBcorr %>%
  mutate(
    Totalbacteria.filter_LB_FBcorr = Totalbacteria.filter_LBcorr - bacteria_FieldBlanks$FieldBlank_TotAvg
  )

bacteria_LBcorr <- left_join(bacteria_LBcorr, SampleInfo, by = c("SampleID", "SampleType"))

bacteria_LBcorr <- bacteria_LBcorr %>%
  mutate(
    RepVolume_L = Volume_L/2, # A and B slide replicates (S9PI vs CW/KOH)
    RepVolume_m3 = RepVolume_L/1000,
    Tot_bacteria.L = Totalbacteria.filter_LBcorr/RepVolume_L,
    Tot_bacteria.m3 = Totalbacteria.filter_LBcorr/RepVolume_m3)

bacteria_ambient <- bacteria_LB_FB_corr %>%
  dplyr::filter(SampleType == "Ambient", Sample != "R6") %>%
  group_by(SampleID) %>%
  summarise(bacteria.m3 = mean(Tot_bacteria.m3),
            sd_bacteria.m3 = sd(Tot_bacteria.m3))

bacteria_bcorr_blue <- bacteria_LB_FB_corr %>%
  dplyr::filter(SampleType != "Ambient" & SampleType != "FieldBlank", Platform == "Blue") %>%
  mutate(
    bcorr_bacteria.m3 = Tot_bacteria.m3 - bacteria_ambient$bacteria.m3,
    Sample_num = as.numeric(str_extract_all(Sample, "\\d+"))) %>%
  select(Sample_num, RepVolume_m3, bcorr_bacteria.m3)



