#########################################################################
###
### Program name: FASMEE23 spore analysis
###
### Purpose: FASMEE23 spore concentrations and EF
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/25/2025
###
#########################################################################

# Decimals in names "." denote "/" 

#Spore counts from Liv Lampman

library(tidyverse)
library(gridExtra)
library(ggstatsplot)


Filter_diam <- 17 ##mm (Obtained from staining tower diameter)
Filter_area <- ((Filter_diam/2)^2)*pi ##mm^2
FOV_diam10x <- 2 ##mm
FOV10x_area <- ((FOV_diam10x/2)^2)*pi ##mm^2
FOV10x.filter <- Filter_area/FOV10x_area

SampleInfo <- read.csv('./FASMEE_Sample_Info20250329.csv', header = T)

SampleInfo <- SampleInfo %>%
  mutate(DATE_UTC = mdy(DATE_UTC, tz = 'MST7MDT'),
         DATE_UTC = ymd(DATE_UTC),
         SampleStart_UTC = ymd_hms(paste(DATE_UTC, SampleStart_UTC), tz = 'UTC'),
         SampleEnd_UTC = ymd_hms(paste(DATE_UTC, SampleEnd_UTC), tz = 'UTC'),
         DateTime_MDT = mdy_hm(DateTime_MDT, tz = 'MST7MDT'),
         DateTime_UTC =  mdy_hm(DateTime_UTC, tz = 'UTC'))

CellSamples <- SampleInfo %>%
  filter(FilterType == 'PTFE', SampleType == 'Smoke' | SampleType == 'Ambient') %>%
  mutate(SampleRep = str_extract(SampleID, ".$"),
         Sample = str_extract(SampleID, "\\d+"))%>%
  filter(SampleRep == "B")

cells <- read.csv('./Cell_Counts_FASMEE23_20250414.csv', header = T)

spores <- cells %>%
  filter(StainType == "CW/KOH") %>%
  select(-LIVECounts, -DEADCounts, -TotalCells) %>%
  rename(SampleID = 'SlideID') %>%
  mutate(SampleID = if_else(SampleType != "LabBlank" & SampleType != "FieldBlank", gsub("_", "", SampleID), SampleID),
         SampleID = if_else(SampleType == "FieldBlank", gsub("_A", "", SampleID), SampleID))

spores <- spores %>%
  mutate(
    LB_Batch = as.factor(case_when(
      StainDate >= 20240503 & StainDate < 20240807 ~ 'A',
      StainDate >= 20240807 & StainDate < 20241119 ~ 'B',
      StainDate >= 20241119 & StainDate < 20241210 ~ 'C',
      ## There were two different TBE solutions with lab blanks from the same stain date, the first set below is TBE1 and second set is TBE2
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank1A_B_20241210", "R5B", "R4B", "R3B", "R2B", "R1B", "B3B", "B2B") ~ 'D',
      StainDate >= 20241210 & StainDate < 20250205 & SampleID %in% 
        c("LabBlank2A_B_20241210", "B6B", "R6B", "R3B", "R7B", "R9B", "B9B", "R10B") ~ 'E', 
      StainDate >= 20250205 ~ 'F',
    )))

spores <- left_join(spores, SampleInfo, by = "SampleID")

spores <- spores %>%
  mutate(
    Volume_L = if_else(is.na(Volume_L), 0, Volume_L),
    RepVolume_L = Volume_L/2, # A and B slide replicates (S9PI vs CW/KOH)
    RepVolume_m3 = RepVolume_L/1000) %>%
  select(-SampleType.y, -Sampler.y) %>%
  rename(SampleType = 'SampleType.x', Sampler = 'Sampler.x')

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
    log1TotalSpores_LB = log1p(TotalSpores_LB),
    TotalSpores_LBcorr = pmax(0, TotalSpores - TotalSpores_LB),
    TotalSpores.filter_LBcorr = TotalSpores_LBcorr*FOV10x.filter,
    TotalSpores_LBcorr_m3 = TotalSpores_LBcorr/RepVolume_m3)

spores <- spores %>%
  mutate(
    Platform = if_else(is.na(Platform), "Blank", Platform),
    SampleType = factor(SampleType),
    Platform = factor(Platform),
    SampleID = factor(SampleID),
    log_volume_offset_m3 = if_else(SampleType == "Smoke" | SampleType == "Ambient", log(RepVolume_m3), 0)) %>% 
  filter(Platform == "Blue" | (Platform == "Red" & SampleType == "Ambient") | Platform == "Blank")

sample_spores <- spores %>%
  group_by(Sample, SampleType, Platform) %>%
  summarise(
    Mean_LB_spores.FOV = mean(TotalSpores_LBcorr),
    Median_LB_spores.FOV = median(TotalSpores_LBcorr),
    Mean_LB_spores.m3 = mean(TotalSpores_LBcorr_m3),
    Median_LB_spores.m3 = median(TotalSpores_LBcorr_m3),
  )

spores_ambient <- spores %>%
  summarise(spores.m3 = mean(Tot_spores.m3),
            sd_spores.m3 = sd(Tot_spores.m3))

spores_bcorr_blue <- spores_LB_FB_corr %>%
  dplyr::filter(SampleType != "Ambient" & SampleType != "FieldBlank", Platform == "Blue") %>%
  mutate(
    bcorr_spores.m3 = Tot_spores.m3 - spores_ambient$spores.m3) %>%
  select(Sample_num, RepVolume_m3, bcorr_spores.m3)  
  
  
  
  
  
  mutate(
    RepVolume_m3 = RepVolume_L/1000,
    Tot_spores.filter = TotalSpores*FOV10x.filter,
    Tot_spores.L = Tot_spores.filter/RepVolume_L,
    Tot_spores.m3 = Tot_spores.filter/RepVolume_m3,
    Project_sample = case_when(Project == "FASMEE" & SampleType == "Smoke" ~ "Forest Smoke",
                               Project == "FASMEE" & SampleType == "Ambient" ~ "Forest Ambient",
                               Project == "Konza4" & SampleType == "Smoke" ~ "Grassland Smoke",
                               Project == "Konza4" & SampleType == "Ambient" ~ "Grassland Ambient"),
    Sample_type = paste0(SampleType, '_', Sample)
  )

k4_spores <- spores %>%
  filter(Project == 'Konza4')

ggplot(spores, aes(Tot_spores.filter)) + geom_boxplot() + 
  facet_wrap(~Sample, scales = 'free') + theme_minimal()

sample_stats <- spores %>%
  group_by(Sample, Project, SampleType, Project_sample) %>%
  summarise(
    meanSpores.filter = mean(Tot_spores.filter),
    medianSpores.filter = median(Tot_spores.filter),
    sdSpores.filter = sd(Tot_spores.filter),
    meanSpores.m3 = mean(Tot_spores.m3),
    medianSpores.m3 = median(Tot_spores.m3),
    sdSpores.m3 = sd(Tot_spores.m3)
  ) %>% ungroup


sample_plt <- sample_stats %>%
  filter(SampleType != "FieldBlank" & SampleType != "LabBlank",
         Project_sample != "Forest Smoke" & Project_sample != "Forest Ambient")

bcorr_sample_plt <- sample_plt %>%
  filter(SampleType != "Ambient") %>%
  mutate(
    bcorr_spores.m3 = case_when(
      Project_sample == "Forest Smoke" ~ medianSpores.m3 - 225,
      Project_sample == "Grassland Smoke" ~ medianSpores.m3 - 18.63),
    bcorr_spores.m3 = if_else(bcorr_spores.m3 < 0, 0, bcorr_spores.m3)
  ) 

Project_plt <- bcorr_sample_plt %>%
  group_by(Project_sample) %>%
  summarise(
    #Spores.filter = median(bcorr_spores.filter),
    Spores.m3 = median(bcorr_spores.m3),
  ) %>% ungroup


plt <- ggplot(sample_plt, aes(x = Project_sample, y = medianSpores.m3)) +
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


plt <- ggbetweenstats(
  data = sample_plt,
  x = Project_sample, 
  y = medianSpores.m3,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

plt <- plt +
  stat_summary(data = Project_plt,
               #aes(label = round(after_stat(y), 2), y = Spores.m3, vjust = -0.5),
               aes(label = round(after_stat(y), 2), y = stage(Spores.m3, after_stat = 175)), 
               fun.y = max, geom = "text", size = 4) + 
  annotate("text", size = 4, x=1, y=175, hjust = 1.5, label = "Median")

plt <- plt + 
  # Add labels and title
  labs(
    x = "",
    y = expression("Spores"~m^-3),
    title = "Wildland Fire Smoke vs. Ambient Air Fungal Spores"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "Roboto", size = 10, color = "black"),
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
plot(plt)

ggsave(
  filename = "./K4_Spore_Boxplt_20241210.png",
  plot = plt,
  width = 7,
  height = 4.5,
  device = "png"
)


