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


######################### Calibration equations
#-----------------------------------------------------------------------------------------

calibration_ints <- read.table('./Input_Data/FASMEE23/FASMEE_Carbon/FASMEE_CO_CO2_Calib_ints.txt', stringsAsFactors = F, sep = ',', header = T)

calibration_ints <- calibration_ints %>% 
  mutate(
    StartTime_MT = ymd_hms(StartTime_MT, tz = 'MST7MDT'),
    EndTime_MT = ymd_hms(EndTime_MT, tz = 'MST7MDT'),
    calib_CO2_ppm = as.numeric(as.character(calib_CO2_ppm)),
    calib_CO_ppm = as.numeric(as.character(calib_CO_ppm)),
    int = seq(1, by = 2, length.out = n()))

calibration_ints_seq <- pivot_longer(calibration_ints, cols = c('StartTime_MT', 'EndTime_MT'), names_to = 'StartStop', values_to = 'Time_MT')

calibration_FASMEE_CO_CO2 <- FASMEE_CO_CO2 %>%
  mutate(int = findInterval(FASMEE_CO_CO2$local_time, calibration_ints_seq$Time_MT)) %>%
  filter(int %% 2 == 1)

calibration_FASMEE_CO_CO2 <- left_join(calibration_FASMEE_CO_CO2, calibration_ints, join_by(int))

# First day calibration 2023 10 08 

calibration_20231008 <- calibration_FASMEE_CO_CO2 %>%
  filter(local_time <= ymd_hms("2023-10-08 23:59:00", tz = "MST7MDT"))

calibration_20231008_CO2 <- calibration_20231008 %>%
  filter(!is.na(calib_CO2_ppm))

calibration_20231008_CO <- calibration_20231008 %>%
  filter(!is.na(calib_CO_ppm))

CO2_20231008 <- lm(calib_CO2_ppm ~ CO2_.ppm., data = calibration_20231008_CO2)
summary(CO2_20231008)$coefficients

d <- expand.grid(CO2_.ppm. = seq(0, 2500), calib_CO2_ppm  = c(408, 1008, 2165))
d$yhat <- predict(CO2_20231008, newdata = d)

p <- ggplot(calibration_20231008_CO2, aes(x = CO2_.ppm., y = calib_CO2_ppm )) + 
  theme_minimal() + geom_point() + geom_line(aes(y = yhat), data = d) 
plot(p)

CO_20231008 <- lm(calib_CO_ppm ~ CO_raw.ppm., data = calibration_20231008_CO)
summary(CO_20231008)$coefficients

d <- expand.grid(CO_raw.ppm. = seq(0, 600), calib_CO_ppm  = c(32, 260, 560))
d$yhat <- predict(CO_20231008, newdata = d)

p <- ggplot(calibration_20231008_CO, aes(x = CO_raw.ppm., y = calib_CO_ppm )) + 
  theme_minimal() + geom_point() + geom_line(aes(y = yhat), data = d) 
plot(p)

# Second day check calibration (2 points) 2023 10 09

calibration_20231009 <- calibration_FASMEE_CO_CO2 %>%
  filter(local_time >= ymd_hms("2023-10-08 23:59:00", tz = "MST7MDT") & local_time <= ymd_hms("2023-10-09 23:59:00", tz = "MST7MDT"))

calibration_20231009_CO2 <- calibration_20231009 %>%
  filter(!is.na(calib_CO2_ppm))

calibration_20231009_CO <- calibration_20231009 %>%
  filter(!is.na(calib_CO_ppm))

CO2_20231009 <- lm(calib_CO2_ppm ~ CO2_.ppm., data = calibration_20231009_CO2)
summary(CO2_20231009)$coefficients

d <- expand.grid(CO2_.ppm. = seq(0, 1500), calib_CO2_ppm  = c(408, 1008))
d$yhat <- predict(CO2_20231009, newdata = d)

p <- ggplot(calibration_20231009_CO2, aes(x = CO2_.ppm., y = calib_CO2_ppm )) + 
  theme_minimal() + geom_point() + geom_line(aes(y = yhat), data = d) 
plot(p)

CO_20231009 <- lm(calib_CO_ppm ~ CO_raw.ppm., data = calibration_20231009_CO)
summary(CO_20231009)$coefficients

d <- expand.grid(CO_raw.ppm. = seq(0, 600), calib_CO_ppm  = c(32, 560))
d$yhat <- predict(CO_20231009, newdata = d)

p <- ggplot(calibration_20231009_CO, aes(x = CO_raw.ppm., y = calib_CO_ppm )) + 
  theme_minimal() + geom_point() + geom_line(aes(y = yhat), data = d) 
plot(p)

# Third day calibration 2023 10 10 

calibration_20231010 <- calibration_FASMEE_CO_CO2 %>%
  filter(local_time >= ymd_hms("2023-10-09 23:59:00", tz = "MST7MDT") & local_time <= ymd_hms("2023-10-11 23:59:00", tz = "MST7MDT"))

calibration_20231010_CO2 <- calibration_20231010 %>%
  filter(!is.na(calib_CO2_ppm))

calibration_20231010_CO <- calibration_20231010 %>%
  filter(!is.na(calib_CO_ppm))

CO2_20231010 <- lm(calib_CO2_ppm ~ CO2_.ppm., data = calibration_20231010_CO2)
summary(CO2_20231010)$coefficients

d <- expand.grid(CO2_.ppm. = seq(0, 1500), calib_CO2_ppm  = c(408, 1008, 2165))
d$yhat <- predict(CO2_20231010, newdata = d)

p <- ggplot(calibration_20231010_CO2, aes(x = CO2_.ppm., y = calib_CO2_ppm )) + 
  theme_minimal() + geom_point() + geom_line(aes(y = yhat), data = d) 
plot(p)

CO_20231010 <- lm(calib_CO_ppm ~ CO_raw.ppm., data = calibration_20231010_CO)
summary(CO_20231010)$coefficients

d <- expand.grid(CO_raw.ppm. = seq(0, 600), calib_CO_ppm  = c(32, 260, 560))
d$yhat <- predict(CO_20231010, newdata = d)

p <- ggplot(calibration_20231010_CO, aes(x = CO_raw.ppm., y = calib_CO_ppm )) + 
  theme_minimal() + geom_point() + geom_line(aes(y = yhat), data = d) 
plot(p)

########################## End of calibration equations