#########################################################################
###
### Program name: FASMEE23 Raw Flight Log Wrangling
###
### Purpose: Adding HAGL and combining F23 flights into one file
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 07/21/2025
###
#########################################################################

### Comments:

# This script uses elevatr::get_elev_point(src = "epqs") which gets elevation data for a given lat/long from a 10 meter DEM

# DJI flight logs also have 10hz readings and they are decreased to 1hz 

# DJI flight logs were uploaded to https://airdata.com for decryption

library(tidyverse)
library(sf)
library(raster)
library(elevatr)
library(geosphere)
library(gt)



################### Blue flight logs

blue_AirData <- list.files('./Input_Data/FASMEE23/FlightLogs/Blue_Lampman/', pattern = 'csv$', full.names = T, recursive = T) %>%
  map(~read.csv(.x)) %>%
  bind_rows()

blue_AirData <- blue_AirData %>%
  mutate(datetime.utc. = ymd_hms(datetime.utc., tz = 'UTC'))

print(nrow(blue_AirData))

blue_AirData <- blue_AirData %>%
  filter(abs(latitude) > 0 & datetime.utc. > ymd('2022-01-01'))

print(nrow(blue_AirData))

blue_AirData <- blue_AirData %>%
  group_by(datetime.utc.) %>% ## Decreasing resolution from 10hz to 1hz
  dplyr::summarise(
    "lat" = mean(latitude),
    "long" = mean(longitude),
    "ASL_m" = mean(altitude_above_seaLevel.feet.)/3.28084,
    "xspeed_m.s" = mean(xSpeed.mph.)/2.237, #not sure which direction is positive and which is negative compared to ignis files
    "yspeed_m.s"= mean(ySpeed.mph.)/2.237,
    "zspeed_m.s" = mean(zSpeed.mph.)/2.237
  ) %>% ungroup

########################## DANGER!!!!!!!!!!!!! #####################
#Running get_elev_point() will download elevations for every lat/long point--lots of time!

# blue_sf <- blue_AirData %>%
#   dplyr::select(
#     lat,
#     long,
#     datetime.utc.
#   )
# 
# blue_sf <- st_as_sf(blue_sf, coords = c('long', 'lat'))
# blue_sf <- st_set_crs(blue_sf, 4326)
# 
# blue_GL<- get_elev_point(blue_sf, prj = 4326, src = "epqs")
# 
# write.csv(blue_GL, './Input_Data/FASMEE23/FlightLogs/blueGL.csv')

##################################################################################################

blue_GL <- read.csv('./Input_Data/FASMEE23/FlightLogs/blueGL.csv', header = F) 

blue_GL = blue_GL[-1, ] %>%
  dplyr::select(V2, V5) %>%
  rename(datetime_utc = "V2", ground_elev_m = "V5") %>%
  mutate(datetime_utc = ymd_hms(datetime_utc, tz = 'UTC'),
         datetime_cdt = with_tz(datetime_utc, tzone = "America/Chicago")
  )

blue_AirData <- left_join(blue_AirData, blue_GL, join_by(datetime.utc. == datetime_utc), copy = TRUE)

blue_AirData <- blue_AirData %>%
  mutate(
    AGL_m = as.numeric(ASL_m) - as.numeric(ground_elev_m)) %>%
  dplyr::select(-ASL_m, -ground_elev_m)

blue_flight_samples <- blue_AirData %>%
  mutate(int = findInterval(blue_AirData$datetime.utc., blue_ints_seq$Time_UTC)) %>%
  filter(int %% 2 == 1)

blue_flight_samples <- left_join(blue_flight_samples, blue_ints, by = "int")

write.csv(blue_flight_samples, './Output/Output_data/FASMEE23/FASMEE23_FlightsAGL.csv', row.names = F)



























