#########################################################################
###
### Program name: Konza 4 PurpleAir 
###
### Purpose: Konza 4 PA analysis
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/06/2025
###
#########################################################################

library(tidyverse)

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


pa_k4 <- pa_k4 %>%
  mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'),
         time_cdt = with_tz(UTCDateTime, tz = "America/Chicago"),
         ## Quality control of particulate matter using the two identical PMS5003 sensors
         difference = abs(pm2_5_atm - pm2_5_atm_b), 
         RPD = abs((pm2_5_atm - pm2_5_atm_b)/((pm2_5_atm + pm2_5_atm_b)/2)*100),
         RPD_diff = if_else(RPD > 70, 1, 0) + if_else(difference > 5, 1, 0),
         QC = case_when(RPD_diff == 2 ~ "bad", TRUE ~ "good"),
         ## Calculating mixing ratio 
         temp_c = (5/9 * (current_temp_f - 32)),
         temp_c = if_else(temp_c < 17, NA, temp_c), ## DO NOT KEEP THIS OUTSIDE OF THIS SCRIPT ################
         dp_c = (5/9 * (current_dewpoint_f - 32)),
         Vapor_pressure = (6.11*10^((7.5*dp_c)/(237.7+dp_c))),
         sat_vapor = (6.11*10^((7.5*temp_c)/(237.7+temp_c))),
         mixing_ratio = (621.97*(Vapor_pressure/(pressure - Vapor_pressure))),
         saturated_mr = (621.97*(sat_vapor/(pressure - sat_vapor))),
         experiment = "k4") %>%
  ## Averages of particulate matter between two sensors
  ## correcting PM 2.5 with EPA formulas 
  mutate(pm1_0_avg = rowMeans(pa_k4[,c('pm1_0_atm', 'pm1_0_atm_b')], na.rm	= TRUE), 
         pm2_5_avg = rowMeans(pa_k4[,c('pm2_5_atm', 'pm2_5_atm_b')], na.rm = TRUE), 
         pm10_0_avg = rowMeans(pa_k4[,c('pm10_0_atm', 'pm10_0_atm_b')], na.rm = TRUE), 
         pm2_5_corr = as.numeric(case_when(
           pm2_5_avg <= 50 ~ ((0.52 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 50 & pm2_5_avg <= 229 ~ ((0.786 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 229 ~ ((0.69 * pm2_5_avg) + ((8.84*10^-4) * (pm2_5_avg)^2) + 2.97))))


samples_k4 <- pa_k4 %>%
  mutate(int = findInterval(pa_k4$time_cdt, leland_ints_seq$time_cdt)) %>%
  filter(int %% 2 == 1)


samples_k4 <- left_join(samples_k4, leland_ints, join_by(int)) %>%
  mutate(pm2_5_corr = ifelse(QC == "bad", NA_character_, pm2_5_corr),
         pm2_5_avg = ifelse(QC == "bad", NA_character_, pm2_5_avg),
         pm10_0_avg = ifelse(QC == "bad", NA_character_, pm10_0_avg),
         pm1_0_avg = ifelse(QC == "bad", NA_character_, pm1_0_avg))


PA_stats_k4 <- samples_k4 %>%
  group_by(Sample) %>%
  dplyr::summarise(
    "MeanPM2.5_ug.m3"  = mean(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "MedianPM2.5_ug.m3"  = median(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "MaxPM2.5_ug.m3" = max(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "MeanPM1.0ug.m3" = mean(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "MedianPM1.0_ug.m3" = median(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "MaxPM1.0_ug.m3" = max(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "MeanPM10_ug.m3" = mean(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "MedianPM10_ug.m3" = median(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "MaxPM10_ug.m3" = max(as.numeric(unlist(pm10_0_avg)), na.rm = T),
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
    SmokeLevel = case_when(
      MedianPM2.5_ug.m3 < 20 ~ "None",
      MedianPM2.5_ug.m3 >= 20 & MedianPM2.5_ug.m3 < 350 ~ "Low",
      MedianPM2.5_ug.m3 >= 350 & MedianPM2.5_ug.m3 < 700 ~ "Moderate",
      MedianPM2.5_ug.m3 >= 700 ~ "High",
      TRUE ~ NA_character_),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")))





