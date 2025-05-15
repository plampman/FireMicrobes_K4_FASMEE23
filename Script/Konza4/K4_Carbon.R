#########################################################################
###
### Program name: Konza 4 EF Carbon Quantification
###
### Purpose: Konza 4 EPA and UI Carbon Quantification
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/06/2025
###
#########################################################################

library(tidyverse)


#list of csv files that represent each sample--taken from "EPA_Carbon_Data_Request_JA.xlsx"
epa_C_files <- list.files('./Input_Data/Konza4/Carbon/Carbon_csv/EPA_smoke/', pattern="\\.csv$", full.names = T)

all_epa_C = list()

#combine all csv files
for (file in epa_C_files) {
  df <- read.csv(file, stringsAsFactors = FALSE, header = T)
  
  df$source <- basename(file)
  
  all_epa_C[[length(all_epa_C) + 1]] <- df
}

epa_C <- do.call(rbind, all_epa_C)

epa_C <- epa_C %>%
  mutate(DateTime_cdt = ymd_hms(DateTime_cdt, tz = "America/Chicago")) %>%
  arrange(DateTime_cdt)

Temp20240408 <- read.csv('./Input_Data/Konza4/Carbon/Temp20240408.csv', header = T)

Temp20240408 <- Temp20240408 %>%
  mutate(DateTime_cdt = ymd_hms(DateTime_cdt, tz = "America/Chicago"))

epa_C <- left_join(epa_C, Temp20240408, by = "DateTime_cdt")

epa_C <- epa_C %>%
  mutate(
    Temp_C = case_when(
      DateTime_cdt < ymd_hms("2024-04-08 23:59:00", tz = "America/Chicago") ~ Temp_C,
      DateTime_cdt > ymd_hms("2024-04-08 23:59:00", tz = "America/Chicago") & 
        DateTime_cdt < ymd_hms("2024-04-09 23:59:00", tz = "America/Chicago") ~ 24,
      DateTime_cdt > ymd_hms("2024-04-09 23:59:00", tz = "America/Chicago") & 
        DateTime_cdt < ymd_hms("2024-04-10 23:59:00", tz = "America/Chicago") ~ 23),
    deltaCO2_ppm = (deltaCO2_mg.m3*(0.082*(273.15+Temp_C)))/44,
    deltaCO_ppm = (deltaCO_mg.m3*(0.082*(273.15+Temp_C)))/28,
    MCE = deltaCO2_ppm/(deltaCO2_ppm + deltaCO_ppm)
  )

samples_epa_C <- epa_C %>%
  mutate(int = findInterval(epa_C$DateTime_cdt, leland_ints_seq$time_cdt)) %>%
  filter(int %% 2 == 1)

samples_epa_C <- left_join(samples_epa_C, leland_ints, join_by(int))

epa_C_mean <- samples_epa_C %>%
  group_by(Sample, source, SampleType) %>%
  summarise(
    MeanCO2_CO_mg.m3 = mean(Carbon_sampled_mg.m3, na.rm = T),
    Median_deltaCO2_CO_mg.m3 = median(Carbon_sampled_mg.m3, na.rm = T),
    Mean_deltaCO2_mg.m3 = mean(deltaCO2_mg.m3, na.rm = T),
    Median_deltaCO_mg.m3 = median(deltaCO_mg.m3, na.rm = T),
    Mean_deltaCO2_ppm = mean(deltaCO2_ppm, na.rm = T),
    Median_deltaCO2_ppm = median(deltaCO2_mg.m3, na.rm = T),
    Mean_deltaCO_ppm = mean(deltaCO_ppm, na.rm = T),
    Median_deltaCO_ppm = median(deltaCO_mg.m3, na.rm = T),
    MeanMCE = mean(MCE, na.rm = T),
    MedianMCE = median(MCE, na.rm = T)
    )

UI_TC <- read.csv('./Input_Data/Konza4/Carbon/Carbon_csv/OCEC.csv', stringsAsFactors = FALSE, header = T)

filter_area_cm2 <- ((3.7/2)^2) * pi #37mm diameter filter for TC

UI_TC <- UI_TC %>%
  mutate(
    TC_mg.filter = (TC_ug.cm2/1000) * filter_area_cm2
  )

UI_TC <- left_join(UI_TC, AL_ints, by = c('SampleID'='Sample'))

UI_TC <- UI_TC %>%
  mutate(
    TC_mg.m3 = TC_mg.filter/Volume_m3
  )

TC_ambient <- UI_TC %>%
  dplyr::filter(SampleType == "A") %>%
  group_by(SampleType) %>%
  summarise(mean_TC_mg.m3 = mean(TC_mg.m3),
            sd_TC_mg.m3 = sd(TC_mg.m3))

UI_TC_bcorr <- UI_TC %>%
  dplyr::filter(SampleType != "A") %>%
  mutate(bcorr_TC_mg.m3 = TC_mg.m3 - TC_ambient$mean_TC_mg.m3)

  
UI_EPA_C <- left_join(epa_C_mean, UI_TC_bcorr, by = c("Sample"="SampleID", "SampleType"))

UI_EPA_C <- UI_EPA_C %>%
  mutate(
    ALL_carbon_mg.m3 = bcorr_TC_mg.m3 + MeanCO2_CO_mg.m3,
    biomass_mg = ALL_carbon_mg.m3*2,
    biomass_g = biomass_mg/1000,
    biomass_kg = biomass_g/1000,
    TC.Gas_C = bcorr_TC_mg.m3/MeanCO2_CO_mg.m3
  ) %>% ungroup

slim_UI_EPA_C <- UI_EPA_C %>%
  select(Sample, 4:13, 31:36)







