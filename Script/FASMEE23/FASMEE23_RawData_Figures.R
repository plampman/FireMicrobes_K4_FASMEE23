#########################################################################
###
### Program name: FASMEE23 raw data figures
###
### Purpose: FASMEE23 microbial analysis
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/14/2025
###
#########################################################################


library(ggplot2)
library(dplyr)
library(readr)
library(stringr)


spores <- spores %>%
  filter(SampleType != "FieldBlank")

# Create a custom function to extract numeric part from SampleID
extract_number <- function(id) {
  num <- as.numeric(str_extract(id, "\\d+"))
  if(is.na(num)) return(9999)  # For non-numeric IDs, assign a high number
  return(num)
}

ordered_data <- spores %>%
  mutate(SampleType = factor(SampleType, levels = c("Ambient", "Smoke"))) %>%
  group_by(SampleID, SampleType) %>%
  summarize(mean_spores = mean(TotalSpores_LBcorr, na.rm = TRUE)) %>%
  mutate(sort_number = extract_number(SampleID)) %>%
  arrange(SampleType, sort_number)

print(ordered_data)

ordered_sample_ids <- ordered_data$SampleID

sample_type_order <- ordered_data %>%
  select(SampleID, SampleType) %>%
  distinct()

separator_positions <- which(diff(as.integer(sample_type_order$SampleType)) != 0) + 0.5

print(paste("Separator positions:", paste(separator_positions, collapse=", ")))

ggplot(spores, aes(x = factor(SampleID, levels = ordered_sample_ids), 
                   y = TotalSpores_LBcorr_m3)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "blue") +
  geom_vline(xintercept = separator_positions, linetype = "dashed", color = "gray50") +
  labs(
    title = "FASMEE23 Spores by FOV",
    x = "Sample ID",
    y = "Total Spores m3"
  ) +
  annotate("rect", 
           xmin = 0.5, 
           xmax = separator_positions[1], 
           ymin = -Inf, 
           ymax = Inf, 
           alpha = 0.1, 
           fill = "lightblue") +
  annotate("rect", 
           xmin = separator_positions[1], 
           xmax = separator_positions[2], 
           ymin = -Inf, 
           ymax = Inf, 
           alpha = 0.1, 
           fill = "lightgreen") +
  annotate("rect", 
           xmin = separator_positions[2], 
           xmax = length(ordered_sample_ids) + 0.5, 
           ymin = -Inf, 
           ymax = Inf, 
           alpha = 0.1, 
           fill = "lightyellow") +
  annotate("text", 
           x = c(separator_positions[1]/2, 
                 (separator_positions[1] + separator_positions[2])/2,
                 (separator_positions[2] + length(ordered_sample_ids) + 0.5)/2), 
           y = Inf, 
           label = c("FieldBlank", "Ambient", "Smoke"),
           vjust = -0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


#ggsave("FASMEE23.m3_boxplot.png", width = 10, height = 6, dpi = 300)
