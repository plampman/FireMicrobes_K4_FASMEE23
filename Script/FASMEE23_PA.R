#########################################################################
###
### Program name: FASMEE23 PurpleAir 
###
### Purpose: FASMEE23 PA analysis
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/25/2025
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


PA_files <- list.files('./Input_Data/FASMEE23/Blue_PA/', pattern="\\.csv$", full.names = T)

all_PA = list()

#combine all csv files
for (file in PA_files) {
  df <- read.csv(file, stringsAsFactors = FALSE, header = T)
  
  df$source <- basename(file)
  
  all_PA[[length(all_PA) + 1]] <- df
}

pa_FASMEE23 <- do.call(rbind, all_PA)


pa_FASMEE23 <- pa_FASMEE23 %>%
  mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'),
         time_mdt = with_tz(UTCDateTime, tz = "MsT7MDT"),
         ## Quality control of particulate matter using the two identical PMS5003 sensors
         difference = abs(pm2_5_atm - pm2_5_atm_b), 
         RPD = abs((pm2_5_atm - pm2_5_atm_b)/((pm2_5_atm + pm2_5_atm_b)/2)*100),
         RPD_diff = if_else(RPD > 70, 1, 0) + if_else(difference > 5, 1, 0),
         QC = case_when(RPD_diff == 2 ~ "bad", TRUE ~ "good"),
         ## Calculating mixing ratio 
         temp_c = (5/9 * (current_temp_f - 32)),
         dp_c = (5/9 * (current_dewpoint_f - 32)),
         Vapor_pressure = (6.11*10^((7.5*dp_c)/(237.7+dp_c))),
         sat_vapor = (6.11*10^((7.5*temp_c)/(237.7+temp_c))),
         mixing_ratio = (621.97*(Vapor_pressure/(pressure - Vapor_pressure))),
         saturated_mr = (621.97*(sat_vapor/(pressure - sat_vapor))),
         experiment = "FASMEE23") %>%
  ## Averages of particulate matter between two sensors
  ## correcting PM 2.5 with EPA formulas 
  mutate(pm1_0_avg = rowMeans(pa_FASMEE23[,c('pm1_0_atm', 'pm1_0_atm_b')], na.rm	= TRUE), 
         pm2_5_avg = rowMeans(pa_FASMEE23[,c('pm2_5_atm', 'pm2_5_atm_b')], na.rm = TRUE), 
         pm10_0_avg = rowMeans(pa_FASMEE23[,c('pm10_0_atm', 'pm10_0_atm_b')], na.rm = TRUE), 
         pm2_5_corr = as.numeric(case_when(
           pm2_5_avg <= 50 ~ ((0.52 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 50 & pm2_5_avg <= 229 ~ ((0.786 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 229 ~ ((0.69 * pm2_5_avg) + ((8.84*10^-4) * (pm2_5_avg)^2) + 2.97))))


PA_samples_FASMEE23 <- pa_FASMEE23 %>%
  mutate(int = findInterval(pa_FASMEE23$time_mdt, blue_ints_seq$Time_UTC)) %>%
  filter(int %% 2 == 1)


PA_samples_FASMEE23 <- left_join(blue_ints, PA_samples_FASMEE23, join_by(int)) %>%
  mutate(pm2_5_corr = ifelse(QC == "bad", NA_character_, pm2_5_corr),
         pm2_5_avg = ifelse(QC == "bad", NA_character_, pm2_5_avg),
         pm10_0_avg = ifelse(QC == "bad", NA_character_, pm10_0_avg),
         pm1_0_avg = ifelse(QC == "bad", NA_character_, pm1_0_avg))

PA_samples_FASMEE23 <- PA_samples_FASMEE23 %>%
  mutate(
    logPM2.5 = if_else(!is.na(pm2_5_corr), log(as.numeric(unlist(pm2_5_corr))), NA),
    sqrtPM2.5 = if_else(!is.na(pm2_5_corr), sqrt(as.numeric(unlist(pm2_5_corr))), NA)
  )


PA_stats_FASMEE23 <- PA_samples_FASMEE23 %>%
  group_by(Sample) %>%
  mutate(Sample = as.numeric(Sample)) %>%
  dplyr::summarise(
    "MeanLogPM2.5" = mean(logPM2.5, na.rm = T),
    "MeansqrtPM2.5" = mean(sqrtPM2.5, na.rm = T),
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
    "MinMR" = min(as.numeric(unlist(mixing_ratio)), na.rm = T)) %>% ungroup

PA_stats_FASMEE23<- PA_stats_FASMEE23 %>%
  mutate(
    SmokeLevel = case_when(
      MedianPM2.5_ug.m3 < 20 ~ "None",
      MedianPM2.5_ug.m3 >= 20 & MedianPM2.5_ug.m3 < 350 ~ "Low",
      MedianPM2.5_ug.m3 >= 350 & MedianPM2.5_ug.m3 < 700 ~ "Moderate",
      MedianPM2.5_ug.m3 >= 700 ~ "High",
      TRUE ~ NA_character_),
    SmokeLevel = factor(SmokeLevel, levels = c("None", "Low", "Moderate", "High")))


stats_FASMEE23 <- samples_FASMEE23 %>%
  group_by(Sample, SampleStart_UTC, SampleEnd_UTC) %>%
  dplyr::summarise(
    "Number Observations" = sum(!is.na(pm2_5_corr)),
    "PA PM2.5 Corrected Mean (ug/m3)"  = mean(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PA PM2.5 Corrected Median (ug/m3)"  = median(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PA PM2.5 Corrected SD (ug/m3)" = sd(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PA PM2.5 Corrected MAX (ug/m3)" = max(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PA PM2.5 Corrected MIN (ug/m3)" = min(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PA PM1.0 Mean (ug/m3)" = mean(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PA PM1.0 Median (ug/m3)" = median(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PA PM1.0 SD (ug/m3)" = sd(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PA PM1.0 Avg. MAX (ug/m3)" = max(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PA PM1.0 Avg. MIN (ug/m3)" = min(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PA PM10 Mean (ug/m3)" = mean(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PA PM10 Median (ug/m3)" = median(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PA PM10 SD (ug/m3)" = sd(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PA PM10 Avg. MAX (ug/m3)" = max(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PA PM10 Avg. MIN (ug/m3)" = min(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PA TEMP (C) MEAN" = mean(as.numeric(unlist(temp_c)), na.rm = T),
    "PA TEMP (C) Median" = median(as.numeric(unlist(temp_c)), na.rm = T),
    "PA TEMP (C) SD" = sd(as.numeric(unlist(temp_c)), na.rm = T),
    "PA TEMP (C) MAX" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "PA TEMP (C) MIN" = min(as.numeric(unlist(temp_c)), na.rm = T),
    "PA % RH MEAN" = mean(as.numeric(unlist(current_humidity)), na.rm = T),
    "PA % RH Median" = median(as.numeric(unlist(current_humidity)), na.rm = T),
    "PA % RH SD" = sd(as.numeric(unlist(current_humidity)), na.rm = T),
    "PA % RH MIN" = min(as.numeric(unlist(current_humidity)), na.rm = T),
    "PA % RH MAX" = max(as.numeric(unlist(current_humidity)), na.rm = T),
    "PA H2O Mixing Ratio MEAN" = mean(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "PA H2O Mixing Ratio Median" = median(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "PA H2O Mixing Ratio SD" = sd(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "PA H2O Mixing Ratio MAX" = max(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "PA H2O Mixing Ratio MIN" = min(as.numeric(unlist(mixing_ratio)), na.rm = T)
  ) %>% arrange(as.numeric(Sample))

PA_sample_box <- ggplot(samples_k4, aes(x = factor(Sample), as.numeric(unlist(pm2_5_corr)))) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7, na.rm = T) + # outlier.shape = NA) +
  #stat_summary(fun.y=mean,aes(shape = 'Mean'), geom="point", size = 1.5, color="blue", fill="blue") +
  #scale_shape_manual('', values = c("Mean" = 18)) + 
  #scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  theme(legend.position="none") 
#coord_cartesian(ylim = c(0, 1.75)) 

PA_sample_box <- PA_sample_box + 
  # Add labels and title
  labs(
    x = "Sample",
    y = expression('PM2.5 '(ug~m^-3)),
  ) + 
  
  theme(
    # This is the new default font in the plot
    text = element_text(family = "sans", size = 8, color = "black"),
    plot.title = element_text(
      family = "sans", 
      size = 14,
      face = "bold",
      color = "#2a475e"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8)
  )

PA_sample_box <- PA_sample_box  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    panel.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"),
    plot.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF")
  )
PA_sample_box




