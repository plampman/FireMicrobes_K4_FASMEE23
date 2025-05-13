

# Tweedie Distribution Assessment for Bacteria Data
# This script analyzes if TotalCells_LBcorr is suitable for Tweedie distribution
# Based on the model: TotalCells_LBcorr ~ SmokeLevel*MedianMR + offset(log_volume_offset_m3) + (1|LB_Batch:SampleID)

# Load required packages
library(tidyverse)    # For data manipulation
library(statmod)      # For tweedie.profile
library(tweedie)      # For Tweedie distribution functions
library(ggplot2)      # For enhanced plotting
library(gridExtra)    # For arranging multiple plots
library(DHARMa)       # For residual diagnostics in GLMMs
library(glmmTMB)      # For fitting the models


# 1. Initial exploration of TotalCells_LBcorr
# ---------------------------------------

# Summary statistics
cat("Summary of TotalCells_LBcorr:\n")
print(summary(bacteria_blue_pa$TotalCells_LBcorr))

# Check for zeros
zero_count <- sum(bacteria_blue_pa$TotalCells_LBcorr == 0, na.rm = TRUE)
zero_percent <- 100 * zero_count / sum(!is.na(bacteria_blue_pa$TotalCells_LBcorr))
cat("\nZeros:", zero_count, "(", round(zero_percent, 2), "%)\n")

# Distribution plots
p1 <- ggplot(bacteria_blue_pa, aes(x = TotalCells_LBcorr)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  labs(title = "Distribution of TotalCells_LBcorr",
       x = "TotalCells_LBcorr", y = "Count") +
  theme_minimal()

# Log-transformed histogram (excluding zeros)
p2 <- ggplot(bacteria_blue_pa %>% filter(TotalCells_LBcorr > 0), 
             aes(x = log10(TotalCells_LBcorr))) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
  labs(title = "Log10 Distribution of TotalCells_LBcorr",
       x = "log10(TotalCells_LBcorr)", y = "Count") +
  theme_minimal()

# Display plots
grid.arrange(p1, p2, ncol = 2)

# 2. Analyze the Variance-Mean Relationship
# ---------------------------------------

# Function to analyze mean-variance relationship by group
analyze_mean_variance <- function(data, response_var, group_var) {
  # Calculate mean and variance within each group
  mean_var_data <- data %>%
    group_by(!!sym(group_var)) %>%
    summarize(
      mean_val = mean(!!sym(response_var), na.rm = TRUE),
      var_val = var(!!sym(response_var), na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 5, mean_val > 0, var_val > 0)  # Filter out groups with insufficient data
  
  # Plot mean vs variance
  p1 <- ggplot(mean_var_data, aes(x = mean_val, y = var_val)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste("Mean-Variance Relationship for", response_var, "by", group_var),
         x = "Group Mean", y = "Group Variance") +
    theme_minimal()
  
  # Log-log plot to estimate Tweedie power parameter
  p2 <- ggplot(mean_var_data, aes(x = log(mean_val), y = log(var_val))) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste("Log-Log Mean-Variance for", response_var, "by", group_var),
         x = "log(Mean)", y = "log(Variance)") +
    theme_minimal()
  
  # Fit the linear model to get the slope (power parameter)
  log_model <- lm(log(var_val) ~ log(mean_val), data = mean_var_data)
  power_estimate <- coef(log_model)[2]
  cat("\nEstimated power parameter for", response_var, "by", group_var, ":", 
      round(power_estimate, 3), "\n")
  cat("R-squared:", round(summary(log_model)$r.squared, 3), "\n\n")
  
  # Print the plots
  grid.arrange(p1, p2, ncol = 2)
  
  return(list(
    power_estimate = power_estimate,
    r_squared = summary(log_model)$r.squared,
    mean_var_data = mean_var_data,
    log_model = log_model
  ))
}

# Analyze variance-mean relationship by different groupings relevant to your model
# 1. By SmokeLevel
mv_smoke <- analyze_mean_variance(bacteria_blue_pa, "TotalCells_LBcorr", "SmokeLevel")

# 2. Create bins for MedianMR (since it's continuous)
bacteria_blue_pa <- bacteria_blue_pa %>%
  mutate(MedianMR_bin = cut_number(MedianMR, n = 10, na.rm = TRUE))
mv_mr <- analyze_mean_variance(bacteria_blue_pa, "TotalCells_LBcorr", "MedianMR_bin")

# 3. By LB_Batch
mv_batch <- analyze_mean_variance(bacteria_blue_pa, "TotalCells_LBcorr", "LB_Batch")

# 4. By SampleType (since it's in your zero-inflation formula)
mv_type <- analyze_mean_variance(bacteria_blue_pa, "TotalCells_LBcorr", "SampleType")

# 3. Check Tweedie Power Parameter with Formal Methods
# --------------------------------------------------

# Function to find optimal Tweedie power parameter
find_optimal_power <- function(data) {
  # Remove rows with NA in response or key predictors
  clean_data <- data %>% 
    filter(!is.na(TotalCells_LBcorr), 
           !is.na(SmokeLevel), 
           !is.na(MedianMR),
           !is.na(log_volume_offset_m3))
  
  # Create a simplified formula for tweedie.profile
  # (without random effects and zero-inflation as tweedie.profile doesn't support them)
  formula_str <- "TotalCells_LBcorr ~ SmokeLevel*MedianMR + offset(log_volume_offset_m3)"
  
  cat("Using formula:", formula_str, "\n")
  
  # Vector of power parameters to check
  p_vec <- seq(1.1, 1.9, by = 0.1)
  
  # Use tweedie.profile to find optimal p
  tryCatch({
    profile <- tweedie.profile(as.formula(formula_str), 
                               p.vec = p_vec, 
                               data = clean_data, 
                               do.plot = TRUE)
    cat("Optimal power parameter:", round(profile$p.max, 3), "\n")
    return(profile)
  }, error = function(e) {
    cat("Error in tweedie.profile:", e$message, "\n")
    return(NULL)
  })
}

# Find optimal power parameter
tweedie_profile <- find_optimal_power(bacteria_blue_pa)

# 4. Model Fitting and Diagnostics
# ------------------------------

# For comparison, fit models with different tweedie power parameters
fit_tweedie_model <- function(data, p_value) {
  clean_data <- data %>% 
    filter(!is.na(TotalCells_LBcorr), 
           !is.na(SmokeLevel), 
           !is.na(MedianMR),
           !is.na(log_volume_offset_m3),
           !is.na(LB_Batch),
           !is.na(SampleID),
           !is.na(SampleType))
  
  # Convert character variables to factors if needed
  if(!is.factor(clean_data$SmokeLevel)) clean_data$SmokeLevel <- as.factor(clean_data$SmokeLevel)
  if(!is.factor(clean_data$SampleType)) clean_data$SampleType <- as.factor(clean_data$SampleType)
  
  # Create model with specified power parameter
  model <- tryCatch({
    glmmTMB(TotalCells_LBcorr ~ SmokeLevel*MedianMR + offset(log_volume_offset_m3) + 
              (1|LB_Batch:SampleID),
            family = tweedie(link = "log", p = p_value), 
            data = clean_data, 
            ziformula = ~SampleType)
  }, error = function(e) {
    cat("Error fitting model with p =", p_value, ":", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(model)) {
    # Print AIC
    cat("AIC for model with p =", p_value, ":", AIC(model), "\n")
    
    # Create DHARMa residuals
    tryCatch({
      res <- simulateResiduals(model)
      plot(res, main = paste("Residual diagnostics (p =", p_value, ")"))
    }, error = function(e) {
      cat("Error creating DHARMa residuals:", e$message, "\n")
    })
  }
  
  return(model)
}

# Try models with different power parameters
# If we found an optimal p from tweedie.profile
if(!is.null(tweedie_profile)) {
  optimal_p <- tweedie_profile$p.max
  # Test a range around the optimal p
  test_p <- c(1.1, optimal_p, 1.9)
} else {
  # Default test values if profile failed
  test_p <- c(1.1, 1.5, 1.9)
}

# Fit models
models <- list()
for(p in test_p) {
  cat("\n--- Fitting model with p =", p, "---\n")
  models[[paste0("p_", p)]] <- fit_tweedie_model(bacteria_blue_pa, p)
}

# 5. Final Assessment
# -----------------

cat("\n\n============= FINAL ASSESSMENT =============\n")

# 1. Check mean-variance relationship power estimates
mv_powers <- c(
  Smoke = mv_smoke$power_estimate,
  MedianMR = mv_mr$power_estimate,
  Batch = mv_batch$power_estimate,
  SampleType = mv_type$power_estimate
)
cat("Mean-variance power estimates:\n")
print(mv_powers)
cat("Average estimated power:", round(mean(mv_powers, na.rm = TRUE), 3), "\n\n")

# 2. Check optimal power from profile
if(!is.null(tweedie_profile)) {
  cat("Tweedie profile optimal power:", round(tweedie_profile$p.max, 3), "\n\n")
} else {
  cat("Tweedie profile analysis failed.\n\n")
}

# 3. Compare models with different power parameters
if(length(models) > 0) {
  aic_values <- sapply(models, function(m) if(!is.null(m)) AIC(m) else NA)
  cat("AIC values for different power parameters:\n")
  print(aic_values)
  
  # Determine best model
  best_p <- names(which.min(aic_values))
  cat("\nBest model by AIC:", best_p, "with AIC =", min(aic_values, na.rm = TRUE), "\n")
}

cat("\nCONCLUSION:\n")
cat("A Tweedie distribution is suitable for modeling TotalCells_LBcorr if:\n")
cat("1. The power parameter estimates are consistently between 1 and 2\n")
cat("2. The mean-variance relationship shows a good fit on log-log scale\n")
cat("3. The model diagnostics (especially DHARMa residuals) look reasonable\n")
cat("4. The data structure matches Tweedie characteristics (zeros + continuous positive values)\n\n")

# Optional: Create a summary visualization of results
# Plot the AIC values for different power parameters if models were successfully fit
if(length(models) > 0 && sum(!is.na(aic_values)) > 1) {
  aic_df <- data.frame(
    power = as.numeric(gsub("p_", "", names(aic_values))),
    AIC = aic_values
  ) %>% filter(!is.na(AIC))
  
  ggplot(aic_df, aes(x = power, y = AIC)) +
    geom_point(size = 3) +
    geom_line() +
    labs(title = "AIC values for different Tweedie power parameters",
         x = "Power parameter (p)", y = "AIC") +
    theme_minimal()
}