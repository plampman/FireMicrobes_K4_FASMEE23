#########################################################################
###
### Program name: Konza 4 spore analysis
###
### Purpose: Konza 4 spore concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 01/29/2025
###
#########################################################################

library(tidyverse)
library(gridExtra)
library(ggstatsplot)

Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam10x <- 2 ##mm
FOV10x_area <- ((FOV_diam10x/2)^2)*pi ##mm^2
FOV10x.filter <- Filter_area/FOV10x_area

SampleInfo <- read.csv('./K4_SampleInfo.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(Sample = as.character(Sample),
         Sample_num = as.numeric(str_extract(Sample, "\\d+")))

cells <- read.csv('./K4_CELL_COUNTS_20250401.csv', header = T)

spores <- cells %>%
  filter(StainType == "CW/KOH") %>%
  select(-LIVECounts, -DEADCounts, -TotalCells, -Clumps.4ormorecells.) %>%
  mutate(Sample = as.character(Sample),
         Sample = str_trim(Sample),
         Sample_num = if_else(SampleType == "Smoke" | SampleType == "Ambient", 
                              as.numeric(str_extract(Sample, "\\d+")), NA_real_))

spores <- left_join(spores, SampleInfo, by = "Sample_num")

spores_plt <- spores %>% filter(SampleType == "Smoke" | SampleType == "Ambient")
ggplot(spores_plt, aes(x = TotalSpores, fill = SampleType)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black") +
  labs(title = "Histogram of Total Spores per FOV by Sample Type",
       x = "Total Spores/FOV",
       y = "Frequency") +
  scale_fill_manual(values = c("Ambient" = "darkgreen", "Smoke" = "orange")) +
  theme_minimal()

spores <- spores %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240621 & StainDate < 20240703 ~ 'A',
      (StainDate >= 20240703 & StainDate < 20241119 & SlideID %in% c("LabBlankE", "12A", "13A", "14A")) ~ 'B',
      (StainDate >= 20240703 & StainDate < 20241119 & !(SlideID %in% c("LabBlankE", "12A", "13A", "14A"))) ~ 'C',
      StainDate >= 20241119 & StainDate < 20250329 ~ 'D',
      StainDate >= 20250329 ~ 'E',
      TRUE ~ NA_character_
    ))
  )

spores <- spores %>%
  mutate(
    RepVolume_m3 = Slide_RepVolume_L/1000
  ) %>%
  select(-Sample.y) %>%
  rename(Sample = 'Sample.x')

blank_means <- spores %>%
  filter(SampleType == "LabBlank") %>%
  group_by(LB_Batch) %>%
  summarise(
    TotalSpores_LB = mean(TotalSpores)
  ) %>%
  ungroup

spores <- left_join(spores, blank_means, by = "LB_Batch")

spores <- spores %>%
  filter(SampleType != "LabBlank") %>%
  mutate(
    TotalSpores_LBcorr = pmax(0, TotalSpores - TotalSpores_LB),
    TotalSpores.filter_LBcorr = TotalSpores_LBcorr*FOV10x.filter,
    TotalSpores_LBcorr_m3 = TotalSpores.filter_LBcorr/RepVolume_m3)

sample_spores <- spores %>%
  group_by(Sample, SampleType) %>%
  summarise(
    Mean_LB_spores.FOV = mean(TotalSpores_LBcorr),
    Median_LB_spores.FOV = median(TotalSpores_LBcorr),
    Mean_LB_spores.m3 = mean(TotalSpores_LBcorr_m3),
    Median_LB_spores.m3 = median(TotalSpores_LBcorr_m3),
  )





plt <- ggplot(sample_plt, aes(x = SampleType, y = Tot_spores.m3)) +
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7, outlier.shape = NA) + 
  geom_jitter(shape = 21, size = 2.5, fill = 'orange', color = 'black', alpha = 0.3, width = 0.15) +
  scale_y_continuous(labels = scales::scientific)

# plt <- plt +
#   stat_summary(data = Project_plt,
#                aes(label = round(after_stat(y), 2), y = Spores.m3, vjust = -0.5),
#                #aes(label = round(after_stat(y), 2), y = stage(Spores.m3, after_stat = 175)), 
#                fun.y = max, geom = "text", size = 4)  
#   #annotate("text", size = 4, x=1, y=175, hjust = 1.5, label = "Median")

plt <- plt + 
  # Add labels and title
  labs(
    x = "",
    y = expression ("Spores "~m^-3),
    title = "Wildland Fire Smoke vs. Ambient Air Fungal Spores"
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
