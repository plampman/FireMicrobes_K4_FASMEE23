#########################################################################
###
### Program name: FASMEE23 qPCR
###
### Purpose: FASMEE23 qPCR wrangling
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 07/16/2025
###
#########################################################################

library(tidyverse)

F23_qPCR <- read_csv('./Input_Data/FASMEE23/FASMEE23_qPCR/FASMEE23_qPCR_20250716.csv') %>%
  left_join(., SampleInfo, by = "SampleID") 

F23_FB_qPCR <- F23_qPCR %>%
  filter(SampleType == "FieldBlank") %>%
  summarise(
    F23_FBmean = mean(DNA18s_copies.sample)
  )

F23_qPCR <- F23_qPCR %>%
  filter(SampleType != "FieldBlank") %>%
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient")) %>%
  mutate(
    Sample_num = as.numeric(str_extract(SampleID, "\\d+")),
    DNA18s_copies_FBcorr.sample = DNA18s_copies.sample - F23_FB_qPCR$F23_FBmean,
    DNA18s_copies_FBcorr.m3 = (DNA18s_copies_FBcorr.sample / (Volume_L/ 1000))
  ) 


write.csv(F23_qPCR, "./Output/Output_data/FASMEE23/FASMEE23_qPCR_SampleInfo.csv", row.names = F)