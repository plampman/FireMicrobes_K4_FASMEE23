## ----- Script to visualize glmmTMB model_logPM_logPM results with interactions and raw data -----
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(scales)
library(patchwork)
library(tidyverse)
library(viridis)

# Function to ensure Sample is properly formatted as a factor
check_and_convert_factors <- function(data) {
  if(!is.factor(data$Sample)) {
    data$Sample <- as.factor(data$Sample)
  }
  return(data)
}

# Apply the function to the data
smoke_spores_pa_C <- check_and_convert_factors(smoke_spores_pa_C)

# Check for NA values and print variable summaries to debug
cat("Summary of logPM25:\n")
print(summary(smoke_spores_pa_C$logPM25))
cat("\nSummary of MedianMCE:\n")
print(summary(smoke_spores_pa_C$MedianMCE))
cat("\nSummary of MedianMR:\n")
print(summary(smoke_spores_pa_C$MedianMR))

# Remove NA values if present
smoke_spores_pa_C <- smoke_spores_pa_C %>% 
  filter(!is.na(logPM25), !is.na(MedianMCE), !is.na(MedianMR))

## ----- Figure 1: Main effect of logPM25 with raw data -----
# Generate predicted values across PM2.5 range using raw data values
log_pm25_range <- range(smoke_spores_pa_C$logPM25, na.rm = TRUE)
log_pm25_seq <- seq(log_pm25_range[1], log_pm25_range[2], length.out = 50) # Reduced points
pm25_seq <- exp(log_pm25_seq)

# Get model_logPM_logPM predictions
em_pm25 <- emmeans(model_logPM_logPM_logPM, ~ logPM25, at = list(logPM25 = log_pm25_seq,
                                                     MedianMCE = mean(smoke_spores_pa_C$MedianMCE, na.rm = TRUE),
                                                     MedianMR = mean(smoke_spores_pa_C$MedianMR, na.rm = TRUE)),
                   offset = 0, # Adjust based on your model_logPM_logPM
                   type = "response")

# Convert to data frame and apply the 72.25 scaling factor
pm25_preds <- as.data.frame(em_pm25)
pm25_preds$PM25 <- pm25_seq
pm25_preds$response <- pm25_preds$response * 72.25
pm25_preds$asymp.LCL <- pm25_preds$asymp.LCL * 72.25
pm25_preds$asymp.UCL <- pm25_preds$asymp.UCL * 72.25

# Create figure with raw data points
fig_pm25 <- ggplot() +
  # Add raw data points
  geom_point(data = smoke_spores_pa_C, 
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores_LBcorr_m3), 
             alpha = 0.4, size = 2, color = "darkgrey") +
  # Add model_logPM_logPM predictions
  geom_ribbon(data = pm25_preds, 
              aes(x = PM25, ymin = asymp.LCL, ymax = asymp.UCL), 
              alpha = 0.2) +
  geom_smooth(data = pm25_preds, 
              aes(x = PM25, y = response),
              method = "loess", se = FALSE, color = "black", linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_x_continuous(trans = "log10") +
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

## ----- Figure 2: Main effect of MedianMCE with raw data -----
# Generate predicted values across MCE range
mce_range <- range(smoke_spores_pa_C$MedianMCE, na.rm = TRUE)
mce_seq <- seq(mce_range[1], mce_range[2], length.out = 50) # Reduced points

# Get model_logPM_logPM predictions
em_mce <- emmeans(model_logPM_logPM_logPM, ~ MedianMCE, at = list(MedianMCE = mce_seq,
                                                      logPM25 = mean(smoke_spores_pa_C$logPM25, na.rm = TRUE),
                                                      MedianMR = mean(smoke_spores_pa_C$MedianMR, na.rm = TRUE)),
                  offset = 0,
                  type = "response")

# Convert to data frame and apply scaling factor
mce_preds <- as.data.frame(em_mce)
mce_preds$response <- mce_preds$response * 72.25
mce_preds$asymp.LCL <- mce_preds$asymp.LCL * 72.25
mce_preds$asymp.UCL <- mce_preds$asymp.UCL * 72.25

# Create figure with raw data points
fig_mce <- ggplot() +
  # Add raw data points
  geom_point(data = smoke_spores_pa_C, 
             aes(x = MedianMCE, y = TotalSpores_LBcorr_m3), 
             alpha = 0.4, size = 2, color = "darkgrey") +
  # Add model_logPM_logPM predictions
  geom_ribbon(data = mce_preds, 
              aes(x = MedianMCE, ymin = asymp.LCL, ymax = asymp.UCL), 
              alpha = 0.2) +
  geom_smooth(data = mce_preds, 
              aes(x = MedianMCE, y = response),
              method = "loess", se = FALSE, color = "black", linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  labs(x = "Modified Combustion Efficiency (MCE)", 
       y = expression(Total~Spores~(spores/m^3)),
       title = "Effect of MCE on Spore Concentrations") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

## ----- Figure 3: Main effect of MedianMR with raw data -----
# Generate predicted values across MR range
mr_range <- range(smoke_spores_pa_C$MedianMR, na.rm = TRUE)
mr_seq <- seq(mr_range[1], mr_range[2], length.out = 50) # Reduced points

# Get model_logPM_logPM predictions
em_mr <- emmeans(model_logPM_logPM_logPM, ~ MedianMR, at = list(MedianMR = mr_seq,
                                                    MedianMCE = mean(smoke_spores_pa_C$MedianMCE, na.rm = TRUE),
                                                    logPM25 = mean(smoke_spores_pa_C$logPM25, na.rm = TRUE)),
                 offset = 0,
                 type = "response")

# Convert to data frame and apply scaling factor
mr_preds <- as.data.frame(em_mr)
mr_preds$response <- mr_preds$response * 72.25
mr_preds$asymp.LCL <- mr_preds$asymp.LCL * 72.25
mr_preds$asymp.UCL <- mr_preds$asymp.UCL * 72.25

# Create figure with raw data points
fig_mr <- ggplot() +
  # Add raw data points
  geom_point(data = smoke_spores_pa_C, 
             aes(x = MedianMR, y = TotalSpores_LBcorr_m3), 
             alpha = 0.4, size = 2, color = "darkgrey") +
  # Add model_logPM_logPM predictions
  geom_ribbon(data = mr_preds, 
              aes(x = MedianMR, ymin = asymp.LCL, ymax = asymp.UCL), 
              alpha = 0.2) +
  geom_smooth(data = mr_preds, 
              aes(x = MedianMR, y = response),
              method = "loess", se = FALSE, color = "black", linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  labs(x = expression(Mixing~Ratio~(g/kg)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = "Effect of Mixing Ratio on Spore Concentrations") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

## ----- Figure 4: Interaction between logPM25 and MedianMCE with raw data -----
# Create grid for different MCE levels - using quantiles to ensure valid values
mce_levels <- quantile(smoke_spores_pa_C$MedianMCE, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

# Get model_logPM_logPM predictions
em_interact_mce <- emmeans(model_logPM_logPM_logPM, ~ logPM25 * MedianMCE, 
                           at = list(logPM25 = log_pm25_seq,
                                     MedianMCE = mce_levels,
                                     MedianMR = mean(smoke_spores_pa_C$MedianMR, na.rm = TRUE)),
                           offset = 0,
                           type = "response")

# Convert to data frame and apply scaling factor
interact_mce_preds <- as.data.frame(em_interact_mce)
interact_mce_preds$PM25 <- rep(pm25_seq, each = length(mce_levels))
interact_mce_preds$response <- interact_mce_preds$response * 72.25
interact_mce_preds$asymp.LCL <- interact_mce_preds$asymp.LCL * 72.25
interact_mce_preds$asymp.UCL <- interact_mce_preds$asymp.UCL * 72.25

# Format MCE values for legend
interact_mce_preds$MCE_formatted <- sprintf("%.2f", interact_mce_preds$MedianMCE)

# Create bins for raw data points to color by MCE proximity
smoke_spores_pa_C <- smoke_spores_pa_C %>%
  mutate(MCE_bin = case_when(
    MedianMCE <= mean(c(mce_levels[1], mce_levels[2])) ~ as.character(round(mce_levels[1], 2)),
    MedianMCE <= mean(c(mce_levels[2], mce_levels[3])) ~ as.character(round(mce_levels[2], 2)),
    TRUE ~ as.character(round(mce_levels[3], 2))
  ))

# Create figure with raw data points
fig_interact_mce <- ggplot() +
  # Add raw data points colored by MCE bin
  geom_point(data = smoke_spores_pa_C, 
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores_LBcorr_m3, color = MCE_bin), 
             alpha = 0.6, size = 2) +
  # Add model_logPM_logPM predictions with smoothing to reduce oscillations
  geom_smooth(data = interact_mce_preds, 
              aes(x = PM25, y = response, color = MCE_formatted, group = MCE_formatted),
              method = "loess", span = 0.5, se = FALSE, linewidth = 1.2) +
  # Add confidence ribbons
  geom_ribbon(data = interact_mce_preds, 
              aes(x = PM25, ymin = asymp.LCL, ymax = asymp.UCL, 
                  fill = MCE_formatted, group = MCE_formatted), 
              alpha = 0.1, color = NA) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_x_continuous(trans = "log10") +
  scale_color_viridis_d(name = "MCE", end = 0.8) +
  scale_fill_viridis_d(name = "MCE", end = 0.8) +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = expression(Interaction~between~PM[2.5]~and~MCE)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA)
  )

## ----- Figure 5: Interaction between logPM25 and MedianMR with raw data -----
# Create grid for different MR levels - using quantiles to ensure valid values
mr_levels <- quantile(smoke_spores_pa_C$MedianMR, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

# Get model_logPM_logPM predictions
em_interact_mr <- emmeans(model_logPM_logPM_logPM, ~ logPM25 * MedianMR, 
                          at = list(logPM25 = log_pm25_seq,
                                    MedianMR = mr_levels,
                                    MedianMCE = mean(smoke_spores_pa_C$MedianMCE, na.rm = TRUE)),
                          offset = 0,
                          type = "response")

# Convert to data frame and apply scaling factor
interact_mr_preds <- as.data.frame(em_interact_mr)
interact_mr_preds$PM25 <- rep(pm25_seq, each = length(mr_levels))
interact_mr_preds$response <- interact_mr_preds$response * 72.25
interact_mr_preds$asymp.LCL <- interact_mr_preds$asymp.LCL * 72.25
interact_mr_preds$asymp.UCL <- interact_mr_preds$asymp.UCL * 72.25

# Format MR values for legend
interact_mr_preds$MR_formatted <- sprintf("%.1f", interact_mr_preds$MedianMR)

# Create bins for raw data points to color by MR proximity
smoke_spores_pa_C <- smoke_spores_pa_C %>%
  mutate(MR_bin = case_when(
    MedianMR <= mean(c(mr_levels[1], mr_levels[2])) ~ as.character(round(mr_levels[1], 1)),
    MedianMR <= mean(c(mr_levels[2], mr_levels[3])) ~ as.character(round(mr_levels[2], 1)),
    TRUE ~ as.character(round(mr_levels[3], 1))
  ))

# Create figure with raw data points
fig_interact_mr <- ggplot() +
  # Add raw data points colored by MR bin
  geom_point(data = smoke_spores_pa_C, 
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores_LBcorr_m3, color = MR_bin), 
             alpha = 0.6, size = 2) +
  # Add model_logPM_logPM predictions with smoothing to reduce oscillations
  geom_smooth(data = interact_mr_preds, 
              aes(x = PM25, y = response, color = MR_formatted, group = MR_formatted),
              method = "loess", span = 0.5, se = FALSE, linewidth = 1.2) +
  # Add confidence ribbons
  geom_ribbon(data = interact_mr_preds, 
              aes(x = PM25, ymin = asymp.LCL, ymax = asymp.UCL, 
                  fill = MR_formatted, group = MR_formatted), 
              alpha = 0.1, color = NA) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_x_continuous(trans = "log10") +
  scale_color_viridis_d(name = "Mixing Ratio", option = "plasma") +
  scale_fill_viridis_d(name = "Mixing Ratio", option = "plasma") +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = expression(Interaction~between~PM[2.5]~and~Mixing~Ratio)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA)
  )

## ----- Combined figure using patchwork -----
# First row: main effects
fig_row1 <- fig_pm25 + fig_mce + fig_mr + 
  plot_layout(ncol = 3, widths = c(1, 1, 1))

# Second row: interaction plots
fig_row2 <- fig_interact_mce + fig_interact_mr + 
  plot_layout(ncol = 2, widths = c(1, 1))

# Combine all figures
combined_fig <- fig_row1 / fig_row2 + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Effects of PM2.5, Combustion Efficiency, and Air Moisture on Spore Concentrations",
    subtitle = "model_logPM_logPM predictions with raw data points",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

# Save figures
ggsave("fig_main_effects_raw.png", fig_row1, width = 15, height = 5, dpi = 300)
ggsave("fig_interactions_raw.png", fig_row2, width = 12, height = 5, dpi = 300)
ggsave("fig_combined_raw.png", combined_fig, width = 15, height = 10, dpi = 300)


### Heat map figure
#--------------------------------------------------------------------------------------------------

## ----- Script to create heatmap visualizing the interaction between logPM25 and MedianMCE -----
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(scales)
library(viridis)
library(tidyverse)

# Increase the emmeans grid size limit
# This allows for more detailed predictions
emm_options(rg.limit = 20000)

# Create a grid of logPM25 and MedianMCE values for heatmap
log_pm25_range <- range(smoke_spores_pa_C$logPM25, na.rm = TRUE)
log_pm25_seq <- seq(log_pm25_range[1], log_pm25_range[2], length.out = 80)
pm25_seq <- exp(log_pm25_seq)

mce_range <- range(smoke_spores_pa_C$MedianMCE, na.rm = TRUE)
mce_seq <- seq(mce_range[1], mce_range[2], length.out = 80)

# Get model_logPM_logPM predictions across the grid
em_grid <- emmeans(model_logPM_logPM_logPM, 
                   ~ logPM25 * MedianMCE, 
                   at = list(
                     logPM25 = log_pm25_seq,
                     MedianMCE = mce_seq,
                     MedianMR = mean(smoke_spores_pa_C$MedianMR, na.rm = TRUE)
                   ),
                   offset = 0,
                   type = "response")

# Convert to data frame
grid_preds <- as.data.frame(em_grid)

# Apply the 72.25 scaling factor
grid_preds$response <- grid_preds$response * 72.25
grid_preds$PM25 <- exp(grid_preds$logPM25)

# Create heatmap
heatmap_plot <- ggplot(grid_preds, aes(x = PM25, y = MedianMCE, fill = response)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    trans = "log10", 
    name = expression(Total~Spores~(spores/m^3)),
    labels = scientific_format(digits = 2),
    option = "plasma"
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(100, 300, 1000, 3000),
    labels = c("100", "300", "1000", "3000")
  ) +
  labs(
    x = expression(PM[2.5]~(μg/m^3)), 
    y = "Modified Combustion Efficiency (MCE)",
    title = expression(Interaction~between~PM[2.5]~and~MCE~on~Spore~Concentrations)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Add contour lines to highlight patterns
heatmap_contour <- heatmap_plot +
  geom_contour(aes(z = response), color = "white", alpha = 0.8, linewidth = 0.5)

# Save the plot
#ggsave("pm25_mce_heatmap.png", heatmap_contour, width = 10, height = 8, dpi = 300)

# Create an enhanced version with raw data points overlaid
heatmap_with_data <- heatmap_contour +
  geom_point(
    data = smoke_spores_pa_C, 
    aes(x = MedianPM2.5_ug.m3, y = MedianMCE),
    size = 2, shape = 21, color = "white", fill = "black", alpha = 0.7
  ) +
  labs(
    subtitle = "White contour lines show equal spore concentration levels, points show raw data"
  )

# Save the enhanced plot
ggsave("pm25_mce_heatmap_with_data.png", heatmap_with_data, width = 10, height = 8, dpi = 300)


# Create a more detailed scatter plot showing each sample uniquely
detailed_scatter <- ggplot(smoke_spores_pa_C, 
                           aes(x = MedianPM2.5_ug.m3, y = TotalSpores_LBcorr_m3, 
                               color = Sample, shape = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "glm", 
              method.args = list(family = "quasipoisson"), 
              color = "black", se = TRUE, linewidth = 1) +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_x_continuous(trans = "log10") +
  scale_color_viridis_d(end = 0.9) +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = "Effect of PM2.5 on Spore Concentrations",
       subtitle = "Each point represents a measurement, color/shape indicate sample origin (n=6)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

# Create faceted plot by sample
facet_by_sample <- ggplot(smoke_spores_pa_C, 
                          aes(x = MedianPM2.5_ug.m3, y = TotalSpores_LBcorr_m3)) +
  geom_point(aes(color = MedianMCE), size = 3) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_x_continuous(trans = "log10") +
  scale_color_viridis_c(name = "MCE") +
  facet_wrap(~ Sample) +
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = "PM2.5 and Spore Concentration by Sample",
       subtitle = "Points colored by Modified Combustion Efficiency") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Combine plots for comparison
combined_plots <- (detailed_scatter / (facet_by_sample + bubble_plot)) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Multiple Ways to Visualize Small Sample Size Data (n=6)",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save the combined plot
ggsave("small_sample_visualization.png", combined_plots, width = 15, height = 12, dpi = 300)


##############################################################################################

library(ggplot2)
library(emmeans)
library(scales)
library(dplyr)

# Assuming your model_logPM is stored as 'model_logPM'

# Step 1: Define a more strategic reference grid with fewer points
# Choose specific PM2.5 values instead of a dense sequence
pm25_values <- seq(min(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                   max(smoke_spores_pa_C$MedianPM2.5_ug.m3), 
                   length.out = 10)  # Reduced number of points

# Choose fewer representative MCE values (e.g., low, medium, high)
mce_values <- quantile(smoke_spores_pa_C$MedianMCE, probs = c(0.1, 0.5, 0.9))

# Fix MedianMR at its median value
median_mr <- mean(smoke_spores_pa_C$MedianMR)

# Create a more manageable prediction grid
pred_grid <- expand.grid(
  logPM25 = log10(pm25_values),
  MedianMCE = mce_values,
  MedianMR = median_mr
)

# Add a typical value for the offset term
pred_grid$log_volume_offset_m3 <- median(smoke_spores_pa_C$log_volume_offset_m3)

# Step 2: Increase the emmeans limit (only if necessary)
emm_options(rg.limit = 30000)  # Increase limit but keep it reasonable

# Step 3: Use emmeans to get predictions
# Specify nuisance parameter to reduce grid size
emm <- emmeans(model_logPM, 
               specs = ~ logPM25 | MedianMCE, 
               at = pred_grid,
               type = "response",
               nuisance = "MedianMR")  # Treating MedianMR as a nuisance parameter

# Convert to a data frame for plotting
pred_df <- as.data.frame(emm)
pred_df$MedianPM2.5_ug.m3 <- 10^pred_df$logPM25

# Now create the plot
interaction_plot <- ggplot() +
  # Add the original data points
  geom_point(data = smoke_spores_pa_C, 
             aes(x = MedianPM2.5_ug.m3, y = TotalSpores_LBcorr_m3, 
                 color = Sample, shape = Sample),
             size = 3, alpha = 0.8) +
  # Add the model_logPM predictions as lines, one for each MCE value
  geom_line(data = pred_df, 
            aes(x = MedianPM2.5_ug.m3, y = response, 
                group = factor(MedianMCE), 
                linetype = factor(round(MedianMCE, 3))),
            color = "black", linewidth = 1) +
  # Add confidence intervals
  geom_ribbon(data = pred_df, 
              aes(x = MedianPM2.5_ug.m3, 
                  ymin = response - SE, 
                  ymax = response + SE,
                  group = factor(MedianMCE)),
              alpha = 0.2) +
  # Use log scales
  scale_y_continuous(trans = "log10", labels = scientific_format(digits = 2)) +
  scale_x_continuous(trans = "log10") +
  # Color and shape for sample points
  scale_color_viridis_d(end = 0.9) +
  # Appropriate labels
  labs(x = expression(PM[2.5]~(μg/m^3)), 
       y = expression(Total~Spores~(spores/m^3)),
       title = "Effect of PM2.5 on Spore Concentrations",
       subtitle = paste("Interaction between PM2.5 and MCE at median MR =", round(median_mr, 3)),
       linetype = "MCE Value") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

print(interaction_plot)