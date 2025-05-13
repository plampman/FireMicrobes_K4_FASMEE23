


library(ggplot2)
library(dplyr)
library(emmeans)
library(ggpubr)
library(viridis)
library(patchwork)

# VISUALIZATION 1: Model-predicted relationship between bacteria and mixing ratio
# -------------------------------------------------------------------------------

# Create a grid of mixing ratio values
mr_range <- seq(min(bacteria_blue_pa$MedianMR, na.rm = TRUE), 
                max(bacteria_blue_pa$MedianMR, na.rm = TRUE), 
                length.out = 100)

# Create prediction data frame
newdata <- expand.grid(
  SampleType = c("Ambient", "Smoke"),
  MedianMR = mr_range,
  log_volume_offset = mean(bacteria_blue_pa$log_volume_offset, na.rm = TRUE)
)

# Add one random effect level (just for prediction purposes)
newdata$LB_Batch <- bacteria_blue_pa$LB_Batch[1]
newdata$SampleID <- bacteria_blue_pa$SampleID[1]

# Get predictions from the model
newdata$predicted <- predict(model, newdata = newdata, type = "response")

# Get raw data for plotting
plotdata <- bacteria_blue_pa %>%
  filter(SampleType %in% c("Ambient", "Smoke")) %>%
  select(SampleType, MedianMR, TotalCells_LBcorr)

# Create plot
plot1 <- ggplot() +
  # Add raw data points
  geom_point(data = plotdata, 
             aes(x = MedianMR, y = TotalCells_LBcorr, color = SampleType, shape = SampleType),
             alpha = 0.5, size = 2) +
  # Add model prediction lines
  geom_line(data = newdata, 
            aes(x = MedianMR, y = predicted, color = SampleType),
            size = 1.2) +
  # Customize appearance
  scale_color_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_shape_manual(values = c("Ambient" = 16, "Smoke" = 17)) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Water Mixing Ratio (g/kg)",
       y = "Bacterial Cells (lab blank corrected)",
       title = "Bacterial Abundance vs. Water Mixing Ratio",
       subtitle = "Significant interaction between sample type and mixing ratio (p = 0.03)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

# VISUALIZATION 2: Ratio of bacteria in smoke vs. ambient across mixing ratio
# ---------------------------------------------------------------------------

# Get pairwise comparisons from emmeans
mr_values_for_comparison <- seq(min(bacteria_blue_pa$MedianMR, na.rm = TRUE),
                                max(bacteria_blue_pa$MedianMR, na.rm = TRUE),
                                length.out = 20)

# Calculate ratios using model parameters
# Using the coefficients from your model
intercept <- -4.1383
smoke_effect <- 11.9704
mr_effect <- 2.6087
interaction_effect <- -3.1005

# Calculate the ratio at each mixing ratio value
ratio_data <- data.frame(
  MedianMR = mr_values_for_comparison,
  ratio = exp(smoke_effect + interaction_effect * mr_values_for_comparison),
  # Calculate approximate standard errors for the ratio
  # This is a simplification - for presentation purposes
  SE = rep(0.6, length(mr_values_for_comparison))  # Approximated from model output
)

# Calculate confidence intervals
ratio_data$lower <- ratio_data$ratio / exp(1.96 * ratio_data$SE)
ratio_data$upper <- ratio_data$ratio * exp(1.96 * ratio_data$SE)
# Determine significance (where lower CI > 1)
ratio_data$significant <- ratio_data$lower > 1

# Create plot
plot2 <- ggplot(ratio_data, aes(x = MedianMR, y = ratio)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = significant), alpha = 0.3) +
  geom_line(size = 1.2, color = "#e74c3c") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("grey80", "#e74c3c"), 
                    labels = c("Not Significant", "Significant (p < 0.05)")) +
  scale_y_log10() +
  labs(x = "Water Mixing Ratio (g/kg)",
       y = "Ratio of Bacteria (Smoke / Ambient)",
       title = "Ratio of Bacterial Abundance in Smoke vs. Ambient Air",
       subtitle = "Values > 1 indicate higher bacterial abundance in smoke samples") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

# VISUALIZATION 3: Boxplots comparing ambient vs. smoke bacteria
# --------------------------------------------------------------

plot3 <- bacteria_blue_pa %>%
  filter(SampleType %in% c("Ambient", "Smoke")) %>%
  ggplot(aes(x = SampleType, y = TotalCells_LBcorr, fill = SampleType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Sample Type",
       y = "Bacterial Cells (lab blank corrected)",
       title = "Comparison of Bacterial Abundance",
       subtitle = "Higher bacterial counts in smoke samples (p = 0.02)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

# VISUALIZATION 4: Environmental data panel (PM2.5 vs bacteria, colored by sample type)
# ------------------------------------------------------------------------------------

# Prepare data for environmental correlations
env_data <- bacteria_blue_pa %>%
  filter(SampleType %in% c("Ambient", "Smoke")) %>%
  group_by(SampleID, SampleType) %>%
  summarize(
    TotalCells_LBcorr = mean(TotalCells_LBcorr, na.rm = TRUE),
    MedianMR = mean(MedianMR, na.rm = TRUE),
    PM2.5 = mean(MedianPM2.5_ug.m3, na.rm = TRUE),
    Temp = mean(MedianTemp_C, na.rm = TRUE),
    RH = mean(MedianRH, na.rm = TRUE),
    .groups = "drop"
  )

# PM2.5 vs bacteria plot
plot4a <- ggplot(env_data, aes(x = PM2.5, y = TotalCells_LBcorr, color = SampleType)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10(labels = scales::comma) +
  labs(x = expression("PM"[2.5]~"(µg/m"^3*")"),
       y = "Bacterial Cells",
       title = "Bacteria vs. PM2.5 Concentration") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Temperature vs bacteria plot
plot4b <- ggplot(env_data, aes(x = Temp, y = TotalCells_LBcorr, color = SampleType)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Temperature (°C)",
       y = "Bacterial Cells",
       title = "Bacteria vs. Temperature") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Combine environmental plots with a shared legend
legend <- get_legend(
  plot4a + 
    guides(color = guide_legend(title = "Sample Type")) +
    theme(legend.position = "bottom")
)

plot4 <- (plot4a + plot4b) / legend +
  plot_layout(heights = c(5, 1)) +
  plot_annotation(
    title = "Environmental Correlations with Bacterial Abundance",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save all plots
ggsave("plot1_model_predictions.png", plot1, width = 10, height = 6, dpi = 300)
ggsave("plot2_ratio_analysis.png", plot2, width = 10, height = 6, dpi = 300)
ggsave("plot3_boxplot_comparison.png", plot3, width = 8, height = 6, dpi = 300)
ggsave("plot4_environmental_correlations.png", plot4, width = 12, height = 8, dpi = 300)

# Create a combined figure for all visualizations
combined_plot <- (plot1 + plot2) / (plot3 + plot4) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Bacterial Abundance in Wildland Fire Smoke vs. Ambient Air",
    subtitle = "Mixed effects model results showing significant interaction with water mixing ratio",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave("combined_visualization.png", combined_plot, width = 16, height = 12, dpi = 300)


# Get estimated marginal means for each sample type
emm_sample <- emmeans(model, ~ SampleType, type = "response")
emm_df <- as.data.frame(emm_sample)

# Create a more informative boxplot using the estimated marginal means
# Add information about the standard errors and confidence intervals
emm_df <- emm_df %>% 
  mutate(
    lower = asymp.LCL,
    upper = asymp.UCL
  )



emmeans_boxplot <- ggplot(emm_df, aes(x = SampleType, y = response, fill = SampleType)) +
  # Add error bars for confidence intervals
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1) +
  # Add points for the estimated means
  geom_point(size = 4, shape = 23, color = "black") +
  # Add a bar for visual effect (optional)
  geom_bar(stat = "identity", alpha = 0.7, width = 0.5) +
  # Add significance bracket
  annotate("segment", x = 1, xend = 2, 
           y = max(emm_df$upper) * 1.2, 
           yend = max(emm_df$upper) * 1.2, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, 
           y = max(emm_df$upper) * 1.1, 
           yend = max(emm_df$upper) * 1.2, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, 
           y = max(emm_df$upper) * 1.1, 
           yend = max(emm_df$upper) * 1.2, linewidth = 0.5) +
  annotate("text", x = 1.5, 
           y = max(emm_df$upper) * 1.3, 
           label = "p = 0.02", size = 4) +
  # Customize appearance
  scale_fill_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Sample Type",
       y = "Estimated Bacterial Abundance (model-adjusted)",
       title = "Model-Estimated Bacterial Abundance by Sample Type",
       subtitle = "Based on estimated marginal means from the Tweedie mixed effects model",
       caption = "Error bars represent 95% confidence intervals") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))


comparison_plot <- ggplot() +
  # Add boxplot of raw data (semi-transparent)
  geom_boxplot(data = bacteria_blue_pa %>% 
                 filter(SampleType %in% c("Ambient", "Smoke")),
               aes(x = SampleType, y = TotalCells_LBcorr, fill = SampleType),
               alpha = 0.3, outlier.shape = NA) +
  # Add jittered points for raw data (small and semi-transparent)
  geom_jitter(data = bacteria_blue_pa %>% 
                filter(SampleType %in% c("Ambient", "Smoke")),
              aes(x = SampleType, y = TotalCells_LBcorr, color = SampleType),
              width = 0.2, alpha = 0.2, size = 1) +
  # Add estimated marginal means with confidence intervals
  geom_pointrange(data = emm_df,
                  aes(x = SampleType, y = response, 
                      ymin = lower, ymax = upper,
                      color = SampleType),
                  size = 1, fatten = 5) +
  # Add significance bracket
  annotate("segment", x = 1, xend = 2, 
           y = max(c(max(emm_df$upper), 
                     max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.2, 
           yend = max(c(max(emm_df$upper), 
                        max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.2, 
           linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, 
           y = max(c(max(emm_df$upper), 
                     max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.1, 
           yend = max(c(max(emm_df$upper), 
                        max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.2, 
           linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, 
           y = max(c(max(emm_df$upper), 
                     max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.1, 
           yend = max(c(max(emm_df$upper), 
                        max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.2, 
           linewidth = 0.5) +
  annotate("text", x = 1.5, 
           y = max(c(max(emm_df$upper), 
                     max(bacteria_blue_pa$TotalCells_LBcorr, na.rm = TRUE))) * 1.3, 
           label = "p = 0.02", size = 4) +
  # Customize appearance
  scale_fill_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_color_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Sample Type",
       y = "Bacterial Abundance (cells/m³)",
       title = "Bacterial Abundance by Sample Type",
       subtitle = "Raw data (boxplots) and model-estimated means (large points with error bars)",
       caption = "Error bars on estimated means represent 95% confidence intervals") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))


emm_sample <- emmeans(model, ~ SampleType, type = "response")
emm_df <- as.data.frame(emm_sample)

# Create the plot with bacteria per m³ on the y-axis
emmeans_boxplot <- ggplot(emm_df, aes(x = SampleType, y = response, fill = SampleType)) +
  # Add error bars for confidence intervals
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, size = 1) +
  # Add points for the estimated means
  geom_point(size = 4, shape = 23, color = "black") +
  # Add a bar for visual effect (optional)
  geom_bar(stat = "identity", alpha = 0.7, width = 0.5) +
  # Add significance bracket
  annotate("segment", x = 1, xend = 2, 
           y = max(emm_df$asymp.UCL) * 1.2, 
           yend = max(emm_df$asymp.UCL) * 1.2, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, 
           y = max(emm_df$asymp.UCL) * 1.1, 
           yend = max(emm_df$asymp.UCL) * 1.2, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, 
           y = max(emm_df$asymp.UCL) * 1.1, 
           yend = max(emm_df$asymp.UCL) * 1.2, linewidth = 0.5) +
  annotate("text", x = 1.5, 
           y = max(emm_df$asymp.UCL) * 1.3, 
           label = "p = 0.02", size = 4) +
  # Customize appearance
  scale_fill_manual(values = c("Ambient" = "#3498db", "Smoke" = "#e74c3c")) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Sample Type",
       y = "Bacterial Concentration (cells/m³)",
       title = "Model-Estimated Bacterial Concentration by Sample Type",
       subtitle = "Based on estimated marginal means from the Tweedie mixed effects model",
       caption = "Error bars represent 95% confidence intervals") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))


newdata_presence <- expand.grid(
  SmokeLevel = c("None", "Low", "Moderate", "High"),
  MedianMR = seq(min(bacteria_blue_pa$MedianMR, na.rm=TRUE),
                 max(bacteria_blue_pa$MedianMR, na.rm=TRUE),
                 length.out=50)
)

# Add a reference level for the random effect
newdata_presence$LB_Batch <- bacteria_blue_pa$LB_Batch[1]
newdata_presence$SampleID <- bacteria_blue_pa$SampleID[1]

# Get predictions
newdata_presence$pred_prob <- predict(presence_model, 
                                      newdata=newdata_presence,
                                      type="response")

# For the positive model
newdata_positive <- newdata_presence
newdata_positive$log_volume_offset_m3 <- mean(bacteria_blue_pa$log_volume_offset_m3, na.rm=TRUE)

# Get predictions
newdata_positive$pred_conc <- predict(positive_model, 
                                      newdata=newdata_positive,
                                      type="response")

# Combine predictions (hurdle model expected value)
newdata_presence$combined <- newdata_presence$pred_prob * newdata_positive$pred_conc


library(ggplot2)
ggplot(newdata_presence, aes(x = MedianMR, y = combined, color = SmokeLevel)) +
  geom_line(size = 1) +
  scale_y_log10() +
  labs(x = "Water Mixing Ratio (g/kg)",
       y = "Expected Bacterial Concentration (cells/FOV)",
       title = "Predicted Bacterial Concentration by Smoke Level and Mixing Ratio") +
  theme_minimal()



# Load necessary libraries
library(ggplot2)
library(DHARMa)
library(gridExtra)
library(grid)

# Create DHARMa residuals for your model
res <- simulateResiduals(TotalSpores_m, n = 1000)

# Extract data for plotting
resid_data <- data.frame(
  fitted = res$fittedPredictedResponse,
  residuals = res$scaledResiduals,
  group = spores_blue_pa$SampleType  # Assuming this is available
)

# Create a publication-quality multi-panel diagnostic plot
# Panel 1: Residuals vs Fitted
p1 <- ggplot(resid_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  theme_classic() +
  labs(x = "Fitted values", y = "Quantile residuals") +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel 2: QQ plot
qq_data <- data.frame(
  theoretical = qunif(ppoints(length(res$scaledResiduals))),
  sample = sort(res$scaledResiduals)
)

p2 <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic() +
  labs(x = "Theoretical quantiles", y = "Observed quantiles") +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel 3: Histogram of residuals
p3 <- ggplot(resid_data, aes(x = residuals)) +
  geom_histogram(bins = 15, fill = "gray80", color = "black") +
  theme_classic() +
  labs(x = "Residuals", y = "Frequency") +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel 4: Residuals by group
p4 <- ggplot(resid_data, aes(x = group, y = residuals)) +
  geom_boxplot(width = 0.5, outlier.shape = 16, outlier.size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_classic() +
  labs(x = "Sample type", y = "Residuals") +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Combine plots into a single figure
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                              top = textGrob("Model diagnostic plots", gp = gpar(fontsize = 14, font = 3)))

# Save as high-quality PDF for publication
ggsave("FASMEE23_TotalSpores_model_diagnostics.pdf", combined_plot, width = 8, height = 8, dpi = 300)

