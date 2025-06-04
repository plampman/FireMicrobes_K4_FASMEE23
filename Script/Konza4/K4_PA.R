#########################################################################
###
### Program name: Konza 4 PurpleAir 
###
### Purpose: Konza 4 PA analysis
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/16/2025
###
#########################################################################

###########################################################################################
### PM2.5 EPA correction equation Sensor Data Cleaning and Correction: 
### Application on the AirNow Fire and Smoke Map 
### 
### In the PDF in the link below on slide 26
### https://cfpub.epa.gov/si/si_public_record_report.cfm?dirEntryId=353088&Lab=CEMM
###########################################################################################

library(tidyverse)

devtools::source_url("https://raw.githubusercontent.com/plampman/PurpleAir_Processing-/main/PurpleAirProcessing.R")

#list of csv files that represent each sample--taken from "EPA_Carbon_Data_Request_JA.xlsx"
PA_files <- list.files('./Input_Data/Konza4/PA/', pattern="\\.csv$", full.names = T)

all_PA = list()

#combine all csv files
for (file in PA_files) {
  df <- read.csv(file, stringsAsFactors = FALSE, header = T)
  
  df$source <- basename(file)
  
  all_PA[[length(all_PA) + 1]] <- df
}

pa_k4 <- do.call(rbind, all_PA)

pa_k4 = process_PA(pa_k4, "America/Chicago", "K4")

samples_k4 <- pa_k4 %>%
  mutate(int = findInterval(pa_k4$local_time, leland_ints_seq$time_cdt)) %>%
  filter(int %% 2 == 1)

samples_k4 <- left_join(samples_k4, leland_ints, join_by(int)) 


PA_stats_k4 <- samples_k4 %>%
  group_by(Sample) %>%
  dplyr::summarise(
    "MeanPM2.5_ug.m3"  = mean(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "MedianPM2.5_ug.m3"  = median(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "MaxPM2.5_ug.m3" = max(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "MeanPM1.0ug.m3" = mean(as.numeric(unlist(pm1_0_atm_avg)), na.rm = T),
    "MedianPM1.0_ug.m3" = median(as.numeric(unlist(pm1_0_atm_avg)), na.rm = T),
    "MaxPM1.0_ug.m3" = max(as.numeric(unlist(pm1_0_atm_avg)), na.rm = T),
    "MeanPM10_ug.m3" = mean(as.numeric(unlist(pm10_0_atm_avg)), na.rm = T),
    "MedianPM10_ug.m3" = median(as.numeric(unlist(pm10_0_atm_avg)), na.rm = T),
    "MaxPM10_ug.m3" = max(as.numeric(unlist(pm10_0_atm_avg)), na.rm = T),
    "MeanTemp_C" = mean(as.numeric(unlist(temp_c)), na.rm = T),
    "MedianTemp_C" = median(as.numeric(unlist(temp_c)), na.rm = T),
    "MaxTemp_C" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "MinTemp_C" = min(as.numeric(unlist(temp_c)), na.rm = T),
    "MeanRH" = mean(as.numeric(unlist(current_humidity)), na.rm = T),
    "MedianRH" = median(as.numeric(unlist(current_humidity)), na.rm = T),
    "MinRH" = min(as.numeric(unlist(current_humidity)), na.rm = T),
    "MaxRH" = max(as.numeric(unlist(current_humidity)), na.rm = T),
    "MeanMR" = mean(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "MedianMR" = median(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "MaxMR" = max(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "MinMR" = min(as.numeric(unlist(mixing_ratio)), na.rm = T)
    ) %>% ungroup


PA_stats_k4 <- PA_stats_k4 %>%
  mutate(
    AQI_PM2.5 = case_when( # AGI levels (using PM concentration) https://www.airnow.gov/aqi/aqi-calculator/
      MedianPM2.5_ug.m3 <= 9 ~ "Good",
      MedianPM2.5_ug.m3 > 9 & MedianPM2.5_ug.m3 <= 125.4 ~ "Moderate/Unhealthy",
      MedianPM2.5_ug.m3 > 125.4 & MedianPM2.5_ug.m3 < 225.4 ~ "Very Unhealthy",
      MedianPM2.5_ug.m3 >= 225.4 ~ "Hazardous",
      TRUE ~ NA_character_),
    AQI_PM2.5 = factor(AQI_PM2.5, 
                       levels = c("Good", "Moderate/Unhealthy", "Very Unhealthy", "Hazardous")),
    SmokeLevel = case_when(
      MedianPM2.5_ug.m3 < 20 ~ "None",
      MedianPM2.5_ug.m3 >= 20 & MedianPM2.5_ug.m3 < 350 ~ "Low",
      MedianPM2.5_ug.m3 >= 350 & MedianPM2.5_ug.m3 < 700 ~ "Moderate",
      MedianPM2.5_ug.m3 >= 700 ~ "High",
      TRUE ~ NA_character_),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")))


# Different "smoke levels"

# Created for our sample pm2.5 distribution (works well with all of our data so far)
#-------------------------------------------------------------------------------------
# SmokeLevel = case_when(
#   MedianPM2.5_ug.m3 < 20 ~ "None",
#   MedianPM2.5_ug.m3 >= 20 & MedianPM2.5_ug.m3 < 350 ~ "Low",
#   MedianPM2.5_ug.m3 >= 350 & MedianPM2.5_ug.m3 < 700 ~ "Moderate",
#   MedianPM2.5_ug.m3 >= 700 ~ "High",
#   TRUE ~ NA_character_),
# SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")))
#------------------------------------------------------------------------------------
# EPA PM2.5 to AQI thresholds
#------------------------------------------------------------------------------------
# AQI_Intervals = case_when(
#   MedianPM2.5_ug.m3 <= 9 ~ "None",
#   MedianPM2.5_ug.m3 > 9 & MedianPM2.5_ug.m3 <= 35.4 ~ "Moderate",
#   MedianPM2.5_ug.m3 > 35.4 & MedianPM2.5_ug.m3 <= 55.4 ~ "UnhealthySensitive",
#   MedianPM2.5_ug.m3 > 55.4 & MedianPM2.5_ug.m3 < 125.4 ~ "Unhealthy",
#   MedianPM2.5_ug.m3 > 125.4 & MedianPM2.5_ug.m3 < 225.4 ~ "VeryUnhealthy",
#   MedianPM2.5_ug.m3 >= 225.4 ~ "Hazardous",
#   TRUE ~ NA_character_),
# AQI_Intervals = factor(AQI_Intervals, 
#                        levels = c("None", "Moderate", "UnhealthySensitive", 
#                                   "Unhealthy", "VeryUnhealthy", "Hazardous")))



    


