#########################################################################
###
### Program name: Konza 4 UAS Locations
###
### Purpose: Konza 4 UAS location comparison EPA and UI
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/06/2025
###
#########################################################################

library(tidyverse)
# require(gt)
library(geosphere)
library(gridExtra)

#Leland pump intervals associated with microbial samples for carbon quantification
leland_ints <- read.csv('./K4_Leland_Intervals.csv', stringsAsFactors = FALSE, header = T)
leland_ints <- leland_ints %>%
  mutate(int = seq(1, by = 2, length.out = n()))

leland_ints_seq <- pivot_longer(leland_ints, 
                                cols = c('Start_DateTime_cdt', 'Stop_DateTime_cdt'), 
                                names_to = 'Start_Stop', values_to = 'time_cdt')
leland_ints_seq <- leland_ints_seq %>%
  mutate(time_cdt = ymd_hms(time_cdt, tz = "America/Chicago"))%>%
  mutate(int = rep(seq(1, by = 2, length.out = ceiling(n()/2)), each = 2))

#list of csv files that represent each sample--taken from "EPA_Carbon_Data_Request_JA.xlsx"
epa_K4_files <- list.files('./Carbon/Carbon_csv/EPA_smoke/', pattern="\\.csv$", full.names = T)

all_epa_K4 = list()

#combine all csv files
for (file in epa_K4_files) {
  df <- read.csv(file, stringsAsFactors = FALSE, header = T)
  
  df$source <- basename(file)
  
  all_epa_K4[[length(all_epa_K4) + 1]] <- df
}

epa_K4 <- do.call(rbind, all_epa_K4)

epa_UAS <- epa_K4 %>%
  mutate(DateTime_cdt = ymd_hms(DateTime_cdt, tz = "America/Chicago")) %>%
  arrange(DateTime_cdt) %>%
  select(DateTime_cdt, Latitude, Longitude, Alt_MSL_ft)

epa_UAS_samples <- epa_UAS %>%
  mutate(int = findInterval(epa_UAS$DateTime_cdt, leland_ints_seq$time_cdt)) %>%
  filter(int %% 2 == 1)

epa_UAS_samples <- left_join(epa_UAS_samples, leland_ints, join_by(int))

epa_UAS_samples <- epa_UAS_samples %>%
  select(DateTime_cdt, Latitude, Longitude, Alt_MSL_ft, Sample, EPA_carbon) %>%
  rename(EPA_Latitude = "Latitude", EPA_Longitude = "Longitude", EPA_MSL_ft = "Alt_MSL_ft")

#UI drone flight records

UI_UAS_files <- list.files('./FlightLogs/AirData_Konza4/', pattern="\\.csv$", full.names = T)

ALL_UAS_files = list()

#combine all csv files
for (file in UI_UAS_files) {
  df <- read.csv(file, stringsAsFactors = FALSE, header = T)
  
  df$source <- basename(file)
  
  ALL_UAS_files[[length(ALL_UAS_files) + 1]] <- df
}

UI_UAS <- do.call(rbind, ALL_UAS_files)

UI_UAS <- UI_UAS %>%
  mutate(datetime.utc. = ymd_hms(datetime.utc., tz = "UTC"),
         DateTime_cdt = with_tz(datetime.utc., tzone = "America/Chicago")) %>%
  arrange(DateTime_cdt)

UI_UAS <- UI_UAS %>%
  group_by(DateTime_cdt) %>%
  summarise(
    UI_Latitude = mean(latitude),
    UI_Longitude = mean(longitude),
    UI_Alt_MSL_ft = mean(altitude_above_seaLevel.feet.)
  )

samples_UI_UAS <- UI_UAS %>%
  mutate(int = findInterval(UI_UAS$DateTime_cdt, leland_ints_seq$time_cdt)) %>%
  filter(int %% 2 == 1)

samples_UI_UAS <- left_join(samples_UI_UAS, leland_ints, join_by(int))

samples_UI_UAS <- samples_UI_UAS %>%
  select(DateTime_cdt, UI_Latitude, UI_Longitude, UI_Alt_MSL_ft, Sample, EPA_carbon)

#Comparison

UI_EPA_UAS <- left_join(samples_UI_UAS, epa_UAS_samples, by = c("DateTime_cdt", "Sample", "EPA_carbon"))
UI_EPA_UAS <- UI_EPA_UAS %>%
  filter(EPA_carbon == TRUE)

UI_EPA_UAS <- UI_EPA_UAS %>%
  mutate(vert_dist_m = abs(UI_Alt_MSL_ft - EPA_MSL_ft)/3.281,
         horizontal_dist_m = distVincentyEllipsoid(
           p1 = cbind(UI_Longitude, UI_Latitude),
           p2 = cbind(EPA_Longitude, EPA_Latitude)),
         angle = bearing(
           p1 = cbind(UI_Longitude, UI_Latitude),
           p2 = cbind(EPA_Longitude, EPA_Latitude)))
         #bearing_by = (angle + 360) %% 360)

UI_EPA_UAS_plt <- UI_EPA_UAS %>%
  filter(Sample == 3)

UI_EPA_UAS_plt <- pivot_longer(
  UI_EPA_UAS_plt, cols = c("vert_dist_m", "horizontal_dist_m"), 
  names_to = "distance_type", values_to = "distance_m")

library(cowplot)

UAS_distance <- ggplot(UI_EPA_UAS_plt, aes(x = DateTime_cdt, y = distance_m, colour = distance_type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("vert_dist_m" = "black", "horizontal_dist_m" = "grey70"),
                     labels = c("UAS Horizontal Distance", "UAS Vertical Distance")) +
  labs(#title = "UIdaho and EPA UAS Distance S26FF",
       #x = "Time CDT",
       y = "UAS Separation (m)",
       color = NULL) + 
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.95),  
    legend.justification = c(1, 1),  
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(6, 6, 6, 6),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = 10)
  )
UAS_distance

plt <- plot_grid(UAS_distance, UI_EPA_carbon, ncol = 1, align = "v", rel_heights = c(1, 1.2))

ggsave('./UIvsEPA_Carbon_S26FF.png',
       plt,
       dpi = 300,
       bg = "white"
       )







stats_yb <- flights_samples_yb %>%
  group_by(sample_comb)%>%
  summarise(
    "Vertical" = mean(vert_dist),
    "Vertical SD" = sd(vert_dist),
    "Horizontal" = mean(horizontal_dist_ft),
    "Horizontal SD" = as.numeric(sd(horizontal_dist_ft)),
    "Blue to Yellow" = mean(bearing_by),
    "Blue to Yellow SD" = sd(bearing_by),
    "Overlap" = max(local_time.x) - min(local_time.x))%>%
  gt(rowname_col = "sample_comb") %>%
  tab_header(
    title = "Blue-Yellow Sample Location Overlap",
    subtitle = "Konza FA and FB November 2021") %>%
  tab_footnote(
    footnote = "Horizontal separation was calculated using Vincenty's Formulae for an ellipsoid",
    locations = cells_body(columns = c("Horizontal", "Horizontal SD"))) %>%
  tab_footnote(
    footnote = "Vertical separation was calculated from the height above ground level, which was derived from the UAS height above sea level and a USGS 1 meter resolution DEM",
    locations = cells_body(columns = c("Vertical", "Vertical SD"))) %>%
  fmt_number(c("Vertical","Vertical SD","Horizontal","Horizontal SD","Blue to Yellow","Blue to Yellow SD"), decimals = 2) %>%
  tab_stubhead(label = "Sample") %>%
  cols_move_to_start("Overlap") %>%
  cols_label(
    "Vertical" = html("Vertical<br>Separation (ft)"),
    "Vertical SD" = html("SD Vertical<br>Separation (ft)"),
    "Horizontal" = html("Horizontal<br>Separation (ft)"),
    "Horizontal SD" = html("SD Horizontal Separation<br>(ft)"),
    "Blue to Yellow" = html("Blue to<br>Yellow (degrees)"),
    "Blue to Yellow SD" = html("SD Blue to<br>Yellow (degrees)"),
    "Overlap" = html("Overlap<br>(seconds)")) %>%
  cols_align(align = "center")

stats_yb




flights_samples_yb <- inner_join(samples_yellow_flights, samples_blue_flights, by = c("local_time" = "local_time"), keep = TRUE)

flights_samples_yb <- flights_samples_yb %>%
  mutate(vert_dist = (height_AGL_ft.x - height_AGL_ft.y),
         horizontal_dist_M = distVincentyEllipsoid(flights_samples_yb[61:60], flights_samples_yb[4:3]),
         horizontal_dist_ft = (horizontal_dist_M*3.281),
         angle = bearing(flights_samples_yb[61:60], flights_samples_yb[4:3]),
         bearing_by = (angle + 360) %% 360)%>%
  unite("sample_comb", sample.x|sample.y, sep = ",", remove = FALSE) 


stats_yb <- flights_samples_yb %>%
  group_by(sample_comb)%>%
  summarise(
    "Vertical" = mean(vert_dist),
    "Vertical SD" = sd(vert_dist),
    "Horizontal" = mean(horizontal_dist_ft),
    "Horizontal SD" = as.numeric(sd(horizontal_dist_ft)),
    "Blue to Yellow" = mean(bearing_by),
    "Blue to Yellow SD" = sd(bearing_by),
    "Overlap" = max(local_time.x) - min(local_time.x))%>%
  gt(rowname_col = "sample_comb") %>%
  tab_header(
    title = "Blue-Yellow Sample Location Overlap",
    subtitle = "Konza FA and FB November 2021") %>%
  tab_footnote(
    footnote = "Horizontal separation was calculated using Vincenty's Formulae for an ellipsoid",
    locations = cells_body(columns = c("Horizontal", "Horizontal SD"))) %>%
  tab_footnote(
    footnote = "Vertical separation was calculated from the height above ground level, which was derived from the UAS height above sea level and a USGS 1 meter resolution DEM",
    locations = cells_body(columns = c("Vertical", "Vertical SD"))) %>%
  fmt_number(c("Vertical","Vertical SD","Horizontal","Horizontal SD","Blue to Yellow","Blue to Yellow SD"), decimals = 2) %>%
  tab_stubhead(label = "Sample") %>%
  cols_move_to_start("Overlap") %>%
  cols_label(
    "Vertical" = html("Vertical<br>Separation (ft)"),
    "Vertical SD" = html("SD Vertical<br>Separation (ft)"),
    "Horizontal" = html("Horizontal<br>Separation (ft)"),
    "Horizontal SD" = html("SD Horizontal Separation<br>(ft)"),
    "Blue to Yellow" = html("Blue to<br>Yellow (degrees)"),
    "Blue to Yellow SD" = html("SD Blue to<br>Yellow (degrees)"),
    "Overlap" = html("Overlap<br>(seconds)")) %>%
  cols_align(align = "center")

stats_yb