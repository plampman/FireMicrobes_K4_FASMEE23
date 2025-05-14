#########################################################################
###
### Program name: FASMEE23 sample and slide metadata
###
### Purpose: FASMEE23 sample and slide metadata
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/13/2025
###
#########################################################################

library(tidyverse)

SampleInfo <- read.csv('./Input_Data/FASMEE23/FASMEE_Sample_Info20250329.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(DATE_UTC = mdy(DATE_UTC, tz = 'MST7MDT'),
         DATE_UTC = ymd(DATE_UTC),
         SampleStart_UTC = ymd_hms(paste(DATE_UTC, SampleStart_UTC), tz = 'UTC'),
         SampleEnd_UTC = ymd_hms(paste(DATE_UTC, SampleEnd_UTC), tz = 'UTC'),
         DateTime_MDT = mdy_hm(DateTime_MDT, tz = 'MST7MDT'),
         DateTime_UTC =  mdy_hm(DateTime_UTC, tz = 'UTC')) 


blue_ints <- SampleInfo %>%
  filter(Platform == 'Blue',
         FilterType == 'PTFE') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+")) %>%
  filter(SampleRep == "A") %>%
  mutate(int = seq(1, by = 2, length.out = n()))

blue_ints_seq <- pivot_longer(blue_ints, cols = c('SampleStart_UTC', 'SampleEnd_UTC'), names_to = 'StartStop', values_to = 'Time_UTC')


CellSamples <- SampleInfo %>%
  filter(FilterType == 'PTFE', SampleType == 'Smoke' | SampleType == 'Ambient') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+"))%>%
  filter(SampleRep == "B")


blue_carbon_ints <- SampleInfo %>%
  filter(Platform == 'Blue',
         Sampler == 'Button' | Sampler == 'CO_CO2') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+")) %>%
  filter(SampleRep == "A") %>%
  mutate(int = seq(1, by = 2, length.out = n()))

blue_quartz <- SampleInfo %>%
  filter(Platform == 'Blue',
         FilterType == 'Quartz') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = as.numeric(str_extract(SampleID, "\\d+")),
         Volume_m3 = Volume_L/1000
  ) 

blue_carbon_ints_seq <- pivot_longer(blue_carbon_ints, cols = c('SampleStart_UTC', 'SampleEnd_UTC'), names_to = 'StartStop', values_to = 'Time_UTC')


SlideInfo <- read.csv('./Input_Data/FASMEE23/FASMEE23_SlideMetadata_20250513.csv', header = T)

