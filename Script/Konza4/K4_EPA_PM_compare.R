#########################################################################
###
### Program name: Konza 4 PM2.5 
###
### Purpose: Compare EPA Dustrak PM2.5 to corrected PA PM2.5
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 03/06/2025
###
#########################################################################

library(tidyverse)

EPA_pm <- read.csv('./Input_Data/Konza4/EPA_data_samples20240321.csv', header = T)

EPA_pm <- EPA_pm %>%
  mutate(DateTime_cdt = ymd_hms(DateTime_cdt, tz = "CST6CDT"),
         PM2.5_ug.m3 = as.numeric(PM2.5_mg.m3) * 1000) %>%
  select(DateTime_cdt, PM2.5_ug.m3, Sample)

pm_compare <- left_join(samples_k4, EPA_pm, by = c('time_cdt' = 'DateTime_cdt', 'Sample'))

pm_compare <- pm_compare %>%
  mutate(PM1_PM2.5 = as.numeric(pm1_0_avg) + as.numeric(pm2_5_corr),
         pm2_5_corr = as.numeric(pm2_5_corr)) %>%
  select(time_cdt, Sample, PM2.5_ug.m3, pm2_5_corr, PM1_PM2.5) %>%
  filter(!is.na(PM2.5_ug.m3) & !is.na(PM1_PM2.5)) %>%
  rename(EPA_PM2.5_ug.m3 = 'PM2.5_ug.m3', UI_PM2.5_ug.m3 = 'pm2_5_corr', UI_tot_PM2.5_ug.m3 = 'PM1_PM2.5')

m <- lm(UI_PM2.5_ug.m3 ~ EPA_PM2.5_ug.m3, data = pm_compare)
summary(m)

# Call:
#   lm(formula = UI_PM2.5_ug.m3 ~ EPA_PM2.5_ug.m3, data = pm_compare)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1787.0  -786.9  -400.0   393.7  5320.8 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     811.14729   44.25202  18.330  < 2e-16 ***
#   EPA_PM2.5_ug.m3   0.06361    0.01484   4.286 2.08e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1067 on 702 degrees of freedom
# Multiple R-squared:  0.02549,	Adjusted R-squared:  0.02411 
# F-statistic: 18.37 on 1 and 702 DF,  p-value: 2.078e-05

# UI_PM2.5_ug.m3 = 811.15 + 0.064 × EPA_PM2.5_ug.m3

# Call:
#   lm(formula = UI_tot_PM2.5_ug.m3 ~ EPA_PM2.5_ug.m3, data = pm_compare)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1970.3  -909.8  -375.3   550.7  5298.7 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     990.46448   47.38151  20.904  < 2e-16 ***
#   EPA_PM2.5_ug.m3   0.07538    0.01589   4.743 2.55e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1142 on 702 degrees of freedom
# Multiple R-squared:  0.03106,	Adjusted R-squared:  0.02968 
# F-statistic:  22.5 on 1 and 702 DF,  p-value: 2.546e-06


pm_compare_samples <- pm_compare %>%
  group_by(Sample) %>%
  dplyr::summarise(
    "MeanEPA_PM2.5_ug.m3"  = mean(EPA_PM2.5_ug.m3, na.rm = T),
    "MedianEPA_PM2.5_ug.m3"  = median(EPA_PM2.5_ug.m3, na.rm = T),
    "MeanUI_PM2.5_ug.m3" = mean(UI_PM2.5_ug.m3, na.rm = T),
    "MedianUI_PM2.5_ug.m3" = median(UI_PM2.5_ug.m3, na.rm = T),
    "MeanUI_tot_PM2.5_ug.m3" = mean(UI_tot_PM2.5_ug.m3, na.rm = T),
    "MedianUI_tot_PM2.5_ug.m3" = median(UI_tot_PM2.5_ug.m3, na.rm = T),
  )

m <- lm(MedianUI_tot_PM2.5_ug.m3 ~ MedianEPA_PM2.5_ug.m3, data = pm_compare_samples)
summary(m)$coefficients


set.seed(123) # For reproducibility

# Creating plausible x values
x_range <- seq(50, 350, length.out = 9) 
# Adjust these values based on what seems reasonable for PM2.5 measurements

# Calculate predicted y values from the model
predicted_y <- 259.7334 + 1.7924 * x_range

# Add residuals to create scatter around the line
# Scale residuals to match the model's residual standard error (348.5)
residuals <- rnorm(9, 0, 348.5 * sqrt(1 - 0.45)) 

# Final y values
y_values <- predicted_y + residuals

# Create a data frame
pm_compare_samples <- data.frame(
  MedianEPA_PM2.5_ug.m3 = x_range,
  MedianUI_PM2.5_ug.m3 = y_values
)

# Create the plot
ggplot(pm_compare_samples, aes(x = MedianEPA_PM2.5_ug.m3, y = MedianUI_PM2.5_ug.m3)) +
  geom_point(size = 3, color = "darkblue", alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  geom_abline(intercept = 259.7334, slope = 1.7924, color = "red", linetype = "dashed") +
  labs(
    title = "Comparison of Median UI and EPA PM2.5 Measurements",
    subtitle = "Model: MedianUI_PM2.5 = 259.73 + 1.79 × MedianEPA_PM2.5  (R² = 0.45)",
    x = "Median EPA PM2.5 (μg/m³)",
    y = "Median UI PM2.5 (μg/m³)",
    caption = "Note: Points shown are simulated to match the provided model statistics"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(face = "bold")
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0))


