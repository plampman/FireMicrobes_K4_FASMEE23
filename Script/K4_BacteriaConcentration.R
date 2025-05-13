#########################################################################
###
### Program name: Konza 4 bacteria analysis
###
### Purpose: Konza 4 bacteria concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/17/2025
###
#########################################################################

library(tidyverse)
library(gridExtra)
library(ggstatsplot)

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam100x <- 0.2 ##mm
FOV100x_area <- ((FOV_diam100x/2)^2)*pi ##mm^2
FOV100x.filter <- Filter_area/FOV100x_area ##FoVs/filter

volume <- read.csv('./Input_Data/Konza4/K4_volume.csv', header = T)

cells <- read.csv('./Input_Data/Konza4/K4_CELL_COUNTS_20250401.csv', header = T)

bacteria <- cells %>%
  filter(StainType == 'S9PI')

sample_bacteria <- bacteria %>%
  group_by(SlideID, Sample, Project, SampleType, StainDate, StainType) %>%
  summarise(
    medianbacteria.FOV = median(TotalCells),
    Tot_bacteria.filter = medianbacteria.FOV*FOV100x.filter,
  ) %>% ungroup

bacteria_LabBlanks <- sample_bacteria %>%
  filter(SampleType == "LabBlank") %>%
  group_by(StainDate) %>%
  summarise(
    "LabBlank_TotAvg" = mean(Tot_bacteria.filter)) %>%
  mutate(StainDate = paste("Blank", StainDate, sep = "_")) %>%
  pivot_wider(names_from = StainDate, values_from = LabBlank_TotAvg)

bacteria_LBcorr <- sample_bacteria %>%
  filter(SampleType != "LabBlank") %>% 
  rowwise() %>%
  mutate(
    Totalbacteria.filter_LBcorr = as.integer(case_when(
      StainDate >= 20240619 & StainDate < 20240703 ~ Tot_bacteria.filter - bacteria_LabBlanks$Blank_20240621,
      StainDate >= 20240703 & StainDate < 20250329 ~ Tot_bacteria.filter - bacteria_LabBlanks$Blank_20240703,
      StainDate >= 20250329 ~ Tot_bacteria.filter - bacteria_LabBlanks$Blank_20250329
    ))) %>%
  ungroup

bacteria_FieldBlanks <- bacteria_LBcorr %>%
  filter(SampleType == "FieldBlank") %>%
  summarise('FieldBlank_TotAvg' = mean(Totalbacteria.filter_LBcorr))

bacteria_FB_LBcorr <- bacteria_LBcorr %>%
  filter(SampleType != 'FieldBlank') %>%
  mutate(Totalbacteria.filter_FB_LBcorr = Totalbacteria.filter_LBcorr - bacteria_FieldBlanks$FieldBlank_TotAvg)

bacteria_FB_LBcorr <- left_join(bacteria_FB_LBcorr, volume, by = "Sample")
  
bacteria_corrected <- bacteria_FB_LBcorr %>%
  mutate(
    RepVolume_m3 = RepVolume_L/1000,
    Tot_bacteria.L = Totalbacteria.filter_FB_LBcorr/RepVolume_L,
    Tot_bacteria.m3 = Totalbacteria.filter_FB_LBcorr/RepVolume_m3)

bacteria_smoke <- bacteria_corrected  %>%
  dplyr::filter(SampleType == "Smoke") %>%
  group_by(SampleType) %>%
  summarise(bacteria.m3 = mean(Tot_bacteria.m3, na.rm = T))

k2_ambient_bacteria.m3 = 285405
k2_ambient_0_bacteria.m3 = 130003
K3_ambient_bacteria.m3 = 338536

bcorr_value = (k2_ambient_0_bacteria.m3 + K3_ambient_bacteria.m3)/2

bacteria_bcorr <- bacteria_corrected %>%
  dplyr::filter(SampleType != "Ambient" & SampleType != "FieldBlank") %>%
  mutate(
    bcorr_bacteria.m3 = Tot_bacteria.m3 - bcorr_value,
    Sample_num = as.numeric(str_extract_all(Sample, "\\d+"))) %>%
  select(Sample_num, RepVolume_m3, Tot_bacteria.m3, bcorr_bacteria.m3)

sample_plt <- bacteria_LBcorr %>%
  filter(SampleType != "FieldBlank")


plt <- ggplot(sample_plt, aes(x = SampleType, y = Tot_bacteria.m3)) +
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7, outlier.shape = NA) + 
  geom_jitter(shape = 21, size = 2.5, fill = 'orange', color = 'black', alpha = 0.3, width = 0.15) +
  scale_y_continuous(labels = scales::scientific)

# plt <- plt +
#   stat_summary(data = Project_plt,
#                aes(label = round(after_stat(y), 2), y = bacteria.m3, vjust = -0.5),
#                #aes(label = round(after_stat(y), 2), y = stage(bacteria.m3, after_stat = 175)), 
#                fun.y = max, geom = "text", size = 4)  
#   #annotate("text", size = 4, x=1, y=175, hjust = 1.5, label = "Median")

plt <- plt + 
  # Add labels and title
  labs(
    x = "",
    y = expression ("bacteria "~m^-3),
    title = "Wildland Fire Smoke vs. Ambient Air Fungal bacteria"
  ) + 
  
  theme(
    # This is the new default font in the plot
    text = element_text(family = "Roboto", size = 8, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 20,
      face = "bold",
      color = "#2a475e"
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 15, 
      face = "bold",
      color="#1b2838"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )
plt
