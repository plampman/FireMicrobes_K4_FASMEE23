#########################################################################
###
### Program name: Konza 4 Microbial emission factors
###
### Purpose: Konza 4 microbial emission quantification
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/06/2025
###
#########################################################################

library(tidyverse)

#Run Konza4_Metadata.R, K4_Carbon.R, and K4_SporeConcentrations.R first

############# Emission factors

UI_EPA_C <- UI_EPA_C %>%
  mutate(Sample = as.numeric(Sample))

Spore_EF <- left_join(UI_EPA_C, spores_bcorr, by = c("Sample"="Sample_num"))

Spore_EF <- Spore_EF %>%
  mutate(spores.kg = bcorr_spores.m3/biomass_kg)

K4_microbe_EF <- left_join(Spore_EF, bacteria_bcorr, by = c("Sample"="Sample_num"))

K4_microbe_EF <- K4_microbe_EF %>%
  mutate(bacteria.kg = bcorr_bacteria.m3/biomass_kg)

K4_Spore_EF <- Spore_EF %>%
  group_by(SampleType) %>%
  summarise(
    q1 = quantile(spores.kg, 0.25),
    medianSporeEF = median(spores.kg),
    q3 = quantile(spores.kg, 0.75),
    iqr = IQR(spores.kg))



ggplot(Spore_EF, aes(x = spores.kg)) +
  geom_histogram(bins = 15, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Spores per Kilogram",
    x = "Spores per Kilogram",
    y = "Frequency"
  ) +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) 