#########################################################################
###
### Program name: Konza4 sample and slide metadata
###
### Purpose: Konza4 sample and slide metadata
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/14/2025
###
#########################################################################


library(tidyverse)

#Leland pump intervals associated with microbial samples for carbon quantification
leland_ints <- read.csv('./Input_Data/Konza4/K4_Leland_Intervals.csv', stringsAsFactors = FALSE, header = T)
leland_ints <- leland_ints %>%
  mutate(int = seq(1, by = 2, length.out = n()))

leland_ints_seq <- pivot_longer(leland_ints, 
                                cols = c('Start_DateTime_cdt', 'Stop_DateTime_cdt'), 
                                names_to = 'Start_Stop', values_to = 'time_cdt')
leland_ints_seq <- leland_ints_seq %>%
  mutate(time_cdt = ymd_hms(time_cdt, tz = "America/Chicago"))%>%
  mutate(int = rep(seq(1, by = 2, length.out = ceiling(n()/2)), each = 2))

#AirLite (total carbon) pump intervals and volume
AL_ints <- read.csv('./Input_Data/Konza4/K4_AirLite_Intervals.csv', stringsAsFactors = FALSE, header = T)

AL_ints <- AL_ints %>%
  mutate(
    Start_DateTime_cdt = ymd_hms(Start_DateTime_cdt, tz = "America/Chicago"),
    Stop_DateTime_cdt = ymd_hms(Stop_DateTime_cdt, tz = "America/Chicago"),
    duration_min = as.numeric(Stop_DateTime_cdt - Start_DateTime_cdt),
    Volume_L = duration_min * AL_FlowRate_L.min,
    Volume_m3 = Volume_L/1000)

SampleInfo <- read.csv('./Input_Data/Konza4/K4_SampleInfo.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(Sample = as.character(Sample),
         Sample_num = as.numeric(str_extract(Sample, "\\d+")))

SlideInfo <- read.csv('./Input_Data/Konza4/Konza4_SlideMetadata_20250513.csv', header = T)


