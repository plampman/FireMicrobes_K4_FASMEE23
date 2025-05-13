## ----- Script to create publication-level figures from glmmTMB model results -----
## Load required packages
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(scales)
library(patchwork)
library(tidyverse)
library(viridis)

# Function to ensure Sample is properly formatted as a factor 
# (in case it wasn't in the original data)
check_and_convert_factors <- function(data) {
  if(!is.factor(data$Sample)) {
    data$Sample <- as.factor(data$Sample)
  }
  return(data)
}

# Apply the function to the data
smoke_spores_pa_C <- check_and_convert_factors(smoke_spores_pa_C)

## ----- Figure 1: Effect of PM2.5 on spores -----
# Generate predicted values across PM2.5 range using emmeans
pm25_seq <- seq(min(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                max(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                length.out = 100)

# Get model predictions
em_pm25 <- emmeans(model_logPM, ~ logPM25, at = list(logPM25 = log(pm25_seq),
                             MedianMCE = mean(smoke_spores_pa_C$MedianMCE),
                             MedianMR = mean(smoke_spores_pa_C$MedianMR)),
                   offset = 0,
                   type = "response")

# Convert to data frame and apply the scaling factor of 72.25
pm25_preds <- as.data.frame(em_pm25)
pm25_preds$MedianPM2.5_ug.m3 <- pm25_seq
pm25_preds$response <- pm25_preds$response * 72.25
pm25_preds$asymp.LCL <- pm25_preds$asymp.LCL * 72.25
pm25_preds$asymp.UCL <- pm25_preds$asymp.UCL * 72.25

# Create figure
fig_pm25 <- ggplot(pm25_preds, aes(x = MedianPM2.5_ug.m3, y = response)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = expression(Effect~of~PM[2.5]~on~Spore~Concentrations)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

## ----- Figure 2: Effect of MCE on spores -----
# Generate predicted values across MCE range
mce_seq <- seq(min(smoke_spores_pa_C$MedianMCE), 
               max(smoke_spores_pa_C$MedianMCE), 
               length.out = 100)

# Get model predictions
em_mce <- emmeans(model_logPM, ~ MedianMCE, at = list(MedianMCE = mce_seq,
                            logPM25 = mean(smoke_spores_pa_C$logPM25),
                            MedianMR = mean(smoke_spores_pa_C$MedianMR)),
                  offset = 0,
                  type = "response")

# Convert to data frame and apply scaling factor
mce_preds <- as.data.frame(em_mce)
mce_preds$response <- mce_preds$response * 72.25
mce_preds$asymp.LCL <- mce_preds$asymp.LCL * 72.25
mce_preds$asymp.UCL <- mce_preds$asymp.UCL * 72.25

# Create figure
fig_mce <- ggplot(mce_preds, aes(x = MedianMCE, y = response)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  labs(x = "Modified Combustion Efficiency (MCE)", 
       y = expression(Total~Spores~(spores/m^3)),
       title = "Effect of MCE on Spore Concentrations") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

## ----- Figure 3: Effect of MR on spores -----
# Generate predicted values across MR range
mr_seq <- seq(min(smoke_spores_pa_C$MedianMR), 
              max(smoke_spores_pa_C$MedianMR), 
              length.out = 100)

# Get model predictions
em_mr <- emmeans(model_logPM, ~ MedianMR, at = list(MedianMR = mr_seq,
                           MedianMCE = mean(smoke_spores_pa_C$MedianMCE),
                           logPM25 = mean(smoke_spores_pa_C$logPM25)),
                 offset = 0,
                 type = "response")

# Convert to data frame and apply scaling factor
mr_preds <- as.data.frame(em_mr)
mr_preds$response <- mr_preds$response * 72.25
mr_preds$asymp.LCL <- mr_preds$asymp.LCL * 72.25
mr_preds$asymp.UCL <- mr_preds$asymp.UCL * 72.25

# Create figure
fig_mr <- ggplot(mr_preds, aes(x = MedianMR, y = response)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  labs(x = expression(Mixing~Ratio~(g/kg)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = "Effect of Mixing Ratio on Spore Concentrations") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

## ----- Figure 4: Interaction between MCE and PM2.5 -----
# Create a grid of MCE and PM2.5 values with wider range of MCE values
# Focus on the low end of MCE as requested
mce_levels <- c(
  min(smoke_spores_pa_C$MedianMCE), 
  min(smoke_spores_pa_C$MedianMCE) + (max(smoke_spores_pa_C$MedianMCE) - min(smoke_spores_pa_C$MedianMCE))/3,
  median(smoke_spores_pa_C$MedianMCE)
)
pm25_seq <- seq(min(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                max(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                length.out = 100)

# Get model predictions for the interaction
em_interact <- emmeans(model_logPM, ~ MedianMCE * logPM25, at = list(MedianMCE = mce_levels,
                                 logPM25 = log(pm25_seq),
                                 MedianMR = mean(smoke_spores_pa_C$MedianMR)),
                       offset = 0,
                       type = "response")

# Convert to data frame and apply scaling factor
interact_preds <- as.data.frame(em_interact)
interact_preds$MedianPM2.5_ug.m3 <- rep(pm25_seq, each = length(mce_levels))
interact_preds$response <- interact_preds$response * 72.25
interact_preds$asymp.LCL <- interact_preds$asymp.LCL * 72.25
interact_preds$asymp.UCL <- interact_preds$asymp.UCL * 72.25

# Format MCE values to show only 2 decimal places for legend
interact_preds$MCE_formatted <- sprintf("%.2f", interact_preds$MedianMCE)

# Create figure
fig_interact <- ggplot(interact_preds, aes(x = MedianPM2.5_ug.m3, 
                                           y = response, 
                                           color = MCE_formatted,
                                           group = MCE_formatted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL, fill = MCE_formatted), 
              alpha = 0.1, color = NA) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_color_viridis_d(name = "MCE", end = 0.8) +
  scale_fill_viridis_d(name = "MCE", end = 0.8) +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = expression(Interaction~between~MCE~and~PM[2.5])) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = c(0.9, 0.85),  # Upper right corner
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA)
  )

## ----- Figure 5: Heatmap of MCE and PM2.5 interaction -----
# Create a grid of MCE and PM2.5 values for heatmap
mce_seq <- seq(min(smoke_spores_pa_C$MedianMCE), 
               max(smoke_spores_pa_C$MedianMCE), 
               length.out = 50)
pm25_seq <- seq(min(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                max(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                length.out = 50)

# Create all combinations
grid_data <- expand.grid(MedianMCE = mce_seq, 
                         MedianPM2.5_ug.m3 = pm25_seq)
grid_data$logPM25 <- log(grid_data$MedianPM2.5_ug.m3)
grid_data$MedianMR <- mean(smoke_spores_pa_C$MedianMR)
grid_data$log_volume_offset_m3 <- 0
# Add Sample variable for random effect (use the first level for predictions)
grid_data$Sample <- levels(smoke_spores_pa_C$Sample)[1]

# Use emmeans instead of direct prediction to handle random effects properly
em_grid <- emmeans(model_logPM, 
                   ~ MedianMCE * logPM25, 
                   at = list(MedianMCE = mce_seq,
                             logPM25 = unique(grid_data$logPM25),
                             MedianMR = mean(smoke_spores_pa_C$MedianMR)),
                   offset = 0,
                   type = "response")

# Convert to data frame
em_grid_df <- as.data.frame(em_grid)

# Match predictions back to grid and apply scaling factor
grid_data <- grid_data %>%
  select(MedianMCE, MedianPM2.5_ug.m3, logPM25) %>%
  left_join(em_grid_df %>% 
              select(MedianMCE, logPM25, response) %>%
              rename(predicted = response),
            by = c("MedianMCE", "logPM25"))

# Apply scaling factor
grid_data$predicted <- grid_data$predicted * 72.25

# Create figure
fig_heatmap <- ggplot(grid_data, aes(x = MedianPM2.5_ug.m3, 
                                     y = MedianMCE, 
                                     fill = predicted)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", 
                       name = expression(Spores/m^3),
                       labels = scientific_format(digits = 2)) +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = "Modified Combustion Efficiency (MCE)",
       title = expression(Spore~Concentration~Response~to~MCE~and~PM[2.5])) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

## ----- Combined figure using patchwork -----
# Keep individual plot legends intact
# No need to collect guides across plots
fig_row1 <- fig_pm25 + fig_mce + fig_mr + 
  plot_layout(ncol = 3)

# Keep individual legends for interaction plots
fig_row2 <- fig_interact + fig_heatmap + 
  plot_layout(ncol = 2)

# Combine all figures
combined_fig <- fig_row1 / fig_row2 + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Effects of PM2.5, Combustion Efficiency, and Air Moisture on Spore Concentrations",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save figures
ggsave("fig_pm25.png", fig_pm25, width = 7, height = 5, dpi = 600)
ggsave("fig_mce.png", fig_mce, width = 7, height = 5, dpi = 600)
ggsave("fig_mr.png", fig_mr, width = 7, height = 5, dpi = 600)
ggsave("fig_interact.png", fig_interact, width = 7, height = 5, dpi = 600)
ggsave("fig_heatmap.png", fig_heatmap, width = 7, height = 5, dpi = 600)
ggsave("K4_SmokeModel.png", combined_fig, width = 15, height = 10, dpi = 600)