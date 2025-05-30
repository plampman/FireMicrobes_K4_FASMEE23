#########################################################################
###
### Program name: FASMEE23 EF Carbon Quantification
###
### Purpose: FASMEE23 EPA and UI Carbon Quantification
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/25/2025
###
#########################################################################

library(tidyverse)


head = read.table('./Input_Data/FASMEE23/FASMEE_Carbon/HEADER.txt', header = FALSE, sep = ",", fill = TRUE) %>%
  filter(!row_number() %in% 1)  
head = as.data.frame(t(head))


#list of csv files that represent each sample--taken from "EPA_Carbon_Data_Request_JA.xlsx"
FASMEE_CO_CO2_files <- list.files('./Input_Data/FASMEE23/FASMEE_Carbon/Data/', pattern="\\.txt$", full.names = T)

all_FASMEE_CO_CO2= list()

#combine all csv files
for (file in FASMEE_CO_CO2_files) {
  df <- read.table(file, stringsAsFactors = FALSE, header = FALSE, sep = ',', fill = TRUE, col.names = head$V1)
  
  df <- df %>% 
    mutate(GPS_Date_Time_.GMT. = mdy_hms(GPS_Date_Time_.GMT., tz = 'UTC'),
           GPS_Date_Time_.GMT. = ymd_hms(GPS_Date_Time_.GMT.)) %>%
    filter(!is.na(GPS_Date_Time_.GMT.))
  
  df$source <- basename(file)
  
  all_FASMEE_CO_CO2[[length(all_FASMEE_CO_CO2) + 1]] <- df
}

FASMEE_CO_CO2 <- do.call(rbind, all_FASMEE_CO_CO2)

FASMEE_CO_CO2 <- FASMEE_CO_CO2 %>%
  mutate(local_time = with_tz(GPS_Date_Time_.GMT., tzone = 'MST7MDT')) %>%
  arrange(local_time)


#Calibrating CO/CO2 data
#---------------------------------------------------------------------------------------

FASMEE_corr_CO_CO2 <- FASMEE_CO_CO2 %>%
  select(source, local_time, GPS_Date_Time_.GMT., CO2_.ppm., CO_raw.ppm., BME_temp.C.)

FASMEE_corr_CO_CO2 <- FASMEE_corr_CO_CO2 %>%
  mutate(
    CO_corr_ppm = case_when(local_time < ymd_hms('2023-10-08 23:59:59', tz = 'MST7MDT') 
                            ~ (1.227681*CO_raw.ppm.) - 9.502209,
                            local_time > ymd_hms('2023-10-08 23:59:59', tz = 'MST7MDT')
                            & local_time < ymd_hms('2023-10-09 23:59:59', tz = 'MST7MDT')
                            ~ (1.203473*CO_raw.ppm.) - 1.339817,
                            local_time > ymd_hms('2023-10-09 23:59:59', tz = 'MST7MDT')
                            & local_time < ymd_hms('2023-10-10 23:59:59', tz = 'MST7MDT')
                            ~ (1.291623*CO_raw.ppm.) - 10.448141),
    CO2_corr_ppm = case_when(local_time < ymd_hms('2023-10-08 23:59:59', tz = 'MST7MDT') 
                            ~ (1.821101*CO2_.ppm.) - 8.289663,
                            local_time > ymd_hms('2023-10-08 23:59:59', tz = 'MST7MDT')
                            & local_time < ymd_hms('2023-10-09 23:59:59', tz = 'MST7MDT')
                            ~ (1.690043*CO2_.ppm.) - 7.265815,
                            local_time > ymd_hms('2023-10-09 23:59:59', tz = 'MST7MDT')
                            & local_time < ymd_hms('2023-10-10 23:59:59', tz = 'MST7MDT')
                            ~ (1.812043*CO2_.ppm.) - 55.479465))

samples_FASMEE_CO_CO2 <- FASMEE_corr_CO_CO2 %>%
  mutate(int = findInterval(FASMEE_CO_CO2$GPS_Date_Time_.GMT., blue_carbon_ints_seq$Time_UTC)) %>%
  filter(int %% 2 == 1)

samples_FASMEE_CO_CO2 <- left_join(samples_FASMEE_CO_CO2, blue_carbon_ints, join_by(int))

unique(samples_FASMEE_CO_CO2$SampleID)

FASMEE_corr_CO_CO2_amb <- samples_FASMEE_CO_CO2 %>%
  filter(SampleType == "Ambient")

ggplot(FASMEE_corr_CO_CO2_amb, aes(x = CO2_corr_ppm)) +
  geom_histogram(bins = 15, fill = "steelblue", color = "black") +
  labs(
    x = "CO2_ppm ambient",
    y = "Frequency") +
  theme_minimal()

amb_thresh_5_6 <- samples_FASMEE_CO_CO2 %>%
  filter(SampleID == "B_pre5_6A") %>%
  summarise(
    CO_mean = mean(CO_corr_ppm),
    CO_sd = sd(CO_corr_ppm),
    CO_thresh = min(CO_corr_ppm),
    CO2_mean = mean(CO2_corr_ppm),
    CO2_median = median(CO2_corr_ppm),
    CO2_sd = sd(CO2_corr_ppm),
    CO2_thresh = CO2_mean + (2*CO2_sd),
    CO2_max = max(CO2_corr_ppm)
  )

amb_thresh_9 <- samples_FASMEE_CO_CO2 %>%
  filter(SampleID == "B_pre9A") %>%
  summarise(
    CO_mean = mean(CO_corr_ppm),
    CO_sd = sd(CO_corr_ppm),
    CO_thresh = min(CO_corr_ppm),
    CO2_mean = mean(CO2_corr_ppm),
    CO2_median = median(CO2_corr_ppm),
    CO2_sd = sd(CO2_corr_ppm),
    CO2_thresh = CO2_mean + (2*CO2_sd),
    CO2_max = max(CO2_corr_ppm)
  )

amb_thresh_10_12 <- samples_FASMEE_CO_CO2 %>%
  filter(SampleID == "B_pre10-12A") %>%
  summarise(
    CO_mean = mean(CO_corr_ppm),
    CO_sd = sd(CO_corr_ppm),
    CO_thresh = min(CO_corr_ppm),
    CO2_mean = mean(CO2_corr_ppm),
    CO2_median = median(CO2_corr_ppm),
    CO2_sd = sd(CO2_corr_ppm),
    CO2_thresh = CO2_mean + (2*CO2_sd),
    CO2_max = max(CO2_corr_ppm)
  )

samples_FASMEE_CO_CO2 <- samples_FASMEE_CO_CO2 %>%
  mutate(
    DeltaCO2_ppm = case_when(
      (SampleID == "B5A" | SampleID == "B6A") & CO2_corr_ppm > amb_thresh_5_6$CO2_thresh ~ CO2_corr_ppm - amb_thresh_5_6$CO2_mean,
      SampleID == "B9A" & CO2_corr_ppm > amb_thresh_9$CO2_thresh ~ CO2_corr_ppm - amb_thresh_9$CO2_mean,
      (SampleID == "B10A" | SampleID == "B11A" | SampleID == "B12A") & CO2_corr_ppm > amb_thresh_10_12$CO2_thresh ~ CO2_corr_ppm - amb_thresh_10_12$CO2_mean,
      TRUE ~ 0),
    DeltaCO2_mg.m3 = (44*DeltaCO2_ppm)/(0.082*(273.15+BME_temp.C.)),
    DeltaCO_ppm = case_when(
      (SampleID == "B5A" | SampleID == "B6A") & DeltaCO2_ppm > 0 ~ CO_corr_ppm - amb_thresh_5_6$CO_mean,
      SampleID == "B9A" & DeltaCO2_ppm > 0 ~ CO_corr_ppm - amb_thresh_9$CO_mean,
      (SampleID == "B10A" | SampleID == "B11A" | SampleID == "B12A") & DeltaCO2_ppm > 0 ~ CO_corr_ppm - amb_thresh_10_12$CO_mean,
      TRUE ~ 0), 
    DeltaCO_mg.m3 = (28*DeltaCO_ppm)/(0.082*(273.15+BME_temp.C.)),
    Gas_Carbon_sampled_mg.m3 = ((12/28)*DeltaCO_mg.m3) + ((12/44)*DeltaCO2_mg.m3),
    MCE = DeltaCO2_ppm/(DeltaCO2_ppm + DeltaCO_ppm)
  )

FASMEE_corr_CO_CO2_mean <- samples_FASMEE_CO_CO2 %>%
  filter(SampleType == "Smoke") %>%
  mutate(Sample = as.numeric(Sample)) %>%
  group_by(Sample, SampleType) %>%
  summarise(
    Mean_CO2_CO_mg.m3 = mean(Gas_Carbon_sampled_mg.m3, na.rm = T),
    Median_CO2_CO_mg.m3 = median(Gas_Carbon_sampled_mg.m3, na.rm = T),
    Mean_DeltaCO2_mg.m3 = mean(DeltaCO2_mg.m3, na.rm = T),
    Median_DeltaCO_mg.m3 = median(DeltaCO_mg.m3, na.rm = T),
    Mean_DeltaCO2_ppm = mean(DeltaCO2_ppm, na.rm = T),
    Median_DeltaCO2_ppm = median(DeltaCO2_mg.m3, na.rm = T),
    Mean_DeltaCO_ppm = mean(DeltaCO_ppm, na.rm = T),
    Median_DeltaCO_ppm = median(DeltaCO_ppm, na.rm = T),
    Max_DeltaCO_ppm = max(DeltaCO_ppm, na.rm = T),
    MeanMCE = mean(MCE, na.rm = T),
    MedianMCE = median(MCE, na.rm = T)
  )

# Total carbon integration
#----------------------------------------------------------------------------------------------------------------

UI_TC <- read.csv('./Input_Data/FASMEE23/FASMEE_Carbon/FASMEE23_TC.csv', stringsAsFactors = FALSE, header = T)

filter_area_cm2 <- ((3.7/2)^2) * pi #37mm diameter filter for TC

UI_TC <- UI_TC %>%
  mutate(
    Sample.ID = as.numeric(str_extract(Sample.ID, "\\d+")),
    TC_mg.filter = (TC_ug.cm2/1000) * filter_area_cm2
  )

UI_TC <- left_join(UI_TC, blue_quartz, by = c('Sample.ID'='Sample'))

UI_TC <- UI_TC %>%
  mutate(
    TC_mg.m3 = TC_mg.filter/Volume_m3
  )

TC_ambient <- UI_TC %>%
  dplyr::filter(SampleType == "Ambient") %>%
  group_by(SampleType) %>%
  summarise(mean_TC_mg.m3 = mean(TC_mg.m3),
            sd_TC_mg.m3 = sd(TC_mg.m3))

UI_TC_bcorr <- UI_TC %>%
  dplyr::filter(SampleType != "Ambient") %>%
  mutate(bcorr_TC_mg.m3 = TC_mg.m3 - TC_ambient$mean_TC_mg.m3)


FASMEE_corr_CO_CO2_mean <- left_join(FASMEE_corr_CO_CO2_mean, UI_TC_bcorr, by = c("Sample"="Sample.ID", "SampleType"))

FASMEE_corr_CO_CO2_mean <- FASMEE_corr_CO_CO2_mean %>%
  mutate(
    Gas_TC_Carbon_mg.m3 = Mean_CO2_CO_mg.m3 + bcorr_TC_mg.m3,
    biomass_mg = Gas_TC_Carbon_mg.m3*2,
    biomass_g = biomass_mg/1000,
    biomass_kg = biomass_g/1000)

slim_fasmmee_C <- FASMEE_corr_CO_CO2_mean %>%
  select(Sample, 3:13, 40:43) %>%
  mutate(
    Sample = as.numeric(Sample)
  )




