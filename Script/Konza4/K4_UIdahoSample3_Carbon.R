#########################################################################
###
### Program name: Konza 4 Uidaho Green box carbon comparison to EPA
###
### Purpose: Konza 4 UI Carbon Quantification
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/06/2025
###
#########################################################################

library(tidyverse)

# Run K4_Carbon.R first


############# UI MPSS box CO/CO2 quantification for 2 ambient samples and 1 smoke sample
# For comparison with EPA CO/CO2 kolibri values

UI_GB_amb <- read.table('./Input_Data/Konza4/GreenBox_20240408_day1/20240408_Ambient1_2.txt', sep = ",", header = T)
UI_GB_smoke <- read.table('./Input_Data/Konza4/GreenBox_20240408_day1/20240408_smoke1.txt', sep = ",", header = T)

UI_GB <- rbind(UI_GB_amb, UI_GB_smoke)

UI_GB <- UI_GB %>%
  mutate(GPS_Date_Time_.GMT. = ymd_hms(GPS_Date_Time_.GMT., tz = 'UTC'),
         local_time = with_tz(GPS_Date_Time_.GMT., tzone = 'CST6CDT'))

samples_UI_GB <- UI_GB %>%
  mutate(int = findInterval(UI_GB$local_time, leland_ints_seq$time_cdt)) %>%
  filter(int %% 2 == 1)

samples_UI_GB <- left_join(samples_UI_GB, leland_ints, join_by(int))

UI_CO_CO2 <- samples_UI_GB %>%
  select(local_time, CO2_.ppm., CO_raw.ppm., Sample, SampleType)

UI_CO_CO2 <- UI_CO_CO2 %>%
  mutate(CO_corr_ppm = (1.3425*CO_raw.ppm.) -6.2514,
         CO2_corr_ppm = (0.9051*CO2_.ppm.) - -253.3502
  )

UI_CO_CO2_amb <- UI_CO_CO2 %>%
  filter(SampleType == "A")

ggplot(UI_CO_CO2_amb, aes(x = CO2_corr_ppm)) +
  geom_histogram(bins = 15, fill = "steelblue", color = "black") +
  labs(
    x = "CO2_ppm ambient",
    y = "Frequency") +
  theme_minimal()

amb_thresh <- UI_CO_CO2 %>%
  filter(SampleType == "A") %>%
  group_by(SampleType) %>%
  summarise(
    CO_mean = mean(CO_corr_ppm),
    CO2_mean = mean(CO2_corr_ppm),
    CO2_sd = sd(CO2_corr_ppm),
    CO2_thresh = CO2_mean + (2*CO2_sd),
    CO2_max = max(CO2_corr_ppm)
  )

PA_samples_K4 <- read.csv('./Input_Data/Konza4/PA_samples_K4.csv', header = T)

PA_samples_K4 <- PA_samples_K4 %>% 
  mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'), 
         time_cdt = ymd_hms(time_cdt, tz = 'CST6CDT'))

PA_join <- PA_samples_K4 %>%
  select(time_cdt, temp_c)

UI_CO_CO2 <- left_join(UI_CO_CO2, PA_join, by = c('local_time' = 'time_cdt'))

library(zoo)

# Simple linear interpolation for PA temperature
UI_CO_CO2$temp_c_filled <- na.approx(UI_CO_CO2$temp_c, na.rm = FALSE)

# Handle any remaining NAs at the beginning or end
UI_CO_CO2$temp_c_filled <- na.locf(UI_CO_CO2$temp_c_filled, fromLast = TRUE, na.rm = FALSE)  # Fill from end
UI_CO_CO2$temp_c_filled <- na.locf(UI_CO_CO2$temp_c_filled, na.rm = FALSE)  # Fill from beginning

UI_CO_CO2$temp_c <- UI_CO_CO2$temp_c_filled
UI_CO_CO2 <- UI_CO_CO2 %>% select(-temp_c_filled)

UI_CO_CO2 <- UI_CO_CO2 %>%
  mutate(
    DeltaCO2_ppm = if_else(CO2_corr_ppm > amb_thresh$CO2_thresh, CO2_corr_ppm - amb_thresh$CO2_mean, 0),
    DeltaCO2_mg.m3 = (44*DeltaCO2_ppm)/(0.082*(273.15+temp_c)),
    CO_mg.m3 = if_else(DeltaCO2_mg.m3 > 0, (28*CO_corr_ppm)/(0.082*(273.15+temp_c)), 0),
    MCE_mg.m3 = DeltaCO2_mg.m3/(DeltaCO2_mg.m3 + CO_mg.m3),
    MCE_ppm = DeltaCO2_ppm/(DeltaCO2_ppm + CO_corr_ppm),
    Carbon_sampled_mg.m3 = ((12/28)*CO_mg.m3) + ((12/44)*DeltaCO2_mg.m3))


UI_CO_CO2_mean <- UI_CO_CO2 %>%
  filter(SampleType == "S") %>%
  group_by(Sample, SampleType) %>%
  summarise(CO2_CO_mean_mg.m3 = mean(Carbon_sampled_mg.m3, na.rm = T),
            biomass_mg = (CO2_CO_mean_mg.m3 + 0.9466879) * 2,
            biomass_g = biomass_mg/1000,
            biomass_kg = biomass_g/1000,
            MCE_mean1 = mean(MCE_mg.m3, na.rm = T),
            MCE_mean2 = mean(MCE_ppm, na.rm = T),
  )

UI_3_EF_spore.kg <- 5562.3329/0.0001761768
epa_3_EF_spore.kg <- 5562.3329/0.0002000954


UI_CO_CO2_plt <- UI_CO_CO2 %>%
  select(local_time, Carbon_sampled_mg.m3, Sample) %>%
  rename(UI_EPA_CO_CO2 = 'Carbon_sampled_mg.m3') %>%
  mutate(source = "UIdaho")

EPA_CO_CO2_plt <- samples_epa_C %>%
  select(DateTime_cdt, Carbon_sampled_mg.m3, Sample) %>%
  rename(UI_EPA_CO_CO2 = 'Carbon_sampled_mg.m3',local_time = 'DateTime_cdt') %>%
  mutate(source = "EPA", UI_EPA_CO_CO2 = replace_na(UI_EPA_CO_CO2, 0))

UI_EPA_plt <- bind_rows(UI_CO_CO2_plt, EPA_CO_CO2_plt)  
UI_EPA_plt <- UI_EPA_plt %>%
  filter(Sample == 3) 

# Create the plot
UI_EPA_carbon <- ggplot(UI_EPA_plt, aes(x = local_time, y = UI_EPA_CO_CO2, color = source)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("UIdaho" = "#1F77B4", "EPA" = "#FF7F0E"),
                     labels = c("EPA Carbon Sampled", "UIdaho Carbon Sampled")) +
  labs(#title = "UIdaho vs. EPA Carbon S26FF",
    x = "Time CDT",
    y = expression("Carbon (mg " * m^{-3} * ")"),
    color = NULL) +  # Setting color = NULL removes the legend title
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.95), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(6, 6, 6, 6),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 10)
  )
UI_EPA_carbon
