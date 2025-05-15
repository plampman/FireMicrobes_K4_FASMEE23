#########################################################################
###
### Program name: Konza 4 model data figures
###
### Purpose: Konza 4 microbial analysis
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 05/15/2025
###
#########################################################################

library(ggplot2)
library(dplyr)
library(emmeans)


### First plot showing spores vs. smoke level - Run K4_Models.R
#---------------------------------------------------------------------------------------

model_predictions <- as.data.frame(summary(em_spores_SmokeLevel, infer = TRUE))

model_predictions <- model_predictions %>%
  rename(
    lowerCI = asymp.LCL,
    upperCI = asymp.UCL
  )

model_predictions$SmokeLevel <- factor(model_predictions$SmokeLevel, 
                                       levels = c("None", "Low", "Moderate", "High"))

raw_data <- spores_pa_C %>%
  group_by(SampleID, SmokeLevel, RepVolume_m3) %>%
  summarize(
    mean_TotalSpores = mean(TotalSpores_LBcorr_m3, na.rm = TRUE),
    .groups = "drop"
  )


p1 <- ggplot() +
  geom_bar(data = model_predictions, 
           aes(x = SmokeLevel, y = response, fill = "Model Prediction"),
           stat = "identity", alpha = 0.6) +
  geom_point(data = model_predictions, 
             aes(x = SmokeLevel, y = response),
             size = 4, shape = 18, color = "black") +
  geom_errorbar(data = model_predictions,
                aes(x = SmokeLevel, ymin = lowerCI, ymax = upperCI),
                width = 0.2, color = "black", size = 1) +
  geom_point(data = raw_data, 
             aes(x = SmokeLevel, y = mean_TotalSpores, color = "Raw Data Means"),
             position = position_jitter(width = 0.2, height = 0), size = 3) +
  labs(title = "Konza Model Predictions and Sample Means",
       subtitle = "Tweedie model with 95% confidence intervals",
       x = "Smoke Level",
       y = expression("Spores m"^-3)) +
  scale_fill_manual(name = "", values = c("Model Prediction" = "#619CFF")) +
  scale_color_manual(name = "", values = c("Raw Data Means" = "#F8766D")) +
  scale_y_continuous(labels = scales::scientific) +
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         color = guide_legend(override.aes = list(size = 3))) +
  theme_bw()

contrast_df <- as.data.frame(summary(contrasts_ratio_spores_SmokeLevel, infer = TRUE))

smoke_levels <- c("None" = 1, "Low" = 2, "Moderate" = 3, "High" = 4)

extract_groups <- function(contrast_name) {
  groups <- strsplit(as.character(contrast_name), " / ")[[1]]
  return(groups)
}

sig_bars <- data.frame(
  comparison = contrast_df$contrast,
  p_value = contrast_df$p.value,
  stringsAsFactors = FALSE
)

sig_bars$groups <- lapply(sig_bars$comparison, extract_groups)
sig_bars$end_group <- sapply(sig_bars$groups, function(x) x[1])  
sig_bars$start_group <- sapply(sig_bars$groups, function(x) x[2]) 

sig_bars$x <- smoke_levels[sig_bars$start_group]
sig_bars$xend <- smoke_levels[sig_bars$end_group]


max_height <- max(c(model_predictions$upperCI, raw_data$mean_TotalSpores), na.rm = TRUE)
sig_bars$y <- max_height * c(1.05, 1.15, 1.25) 

sig_bars$label_text <- vapply(seq_len(nrow(sig_bars)), function(i) {
  p <- sig_bars$p_value[i]
  
  if (p < 0.001) {
    mantissa <- round(p * 10^(-floor(log10(p))), 1)  
    exponent <- floor(log10(p))                      
    
    p_formatted <- paste0("paste(", mantissa, ", \"e\", ", exponent, ")")
  } else if (p < 0.01) {
    p_formatted <- as.character(round(p, 4))
  } else if (p == 0) {
    p_formatted <- "0.00"
  } else {
    p_formatted <- as.character(round(p, 2))
  }
  
  if (p < 0.001) {
    return(paste0("\"***\"~p[Dunnett-adj.] == ", p_formatted))
  } else if (p < 0.01) {
    return(paste0("\"**\"~p[Dunnett-adj.] == ", p_formatted))
  } else if (p < 0.05) {
    return(paste0("\"*\"~p[Dunnett-adj.] == ", p_formatted))
  } else {
    return(paste0("p[Dunnett-adj.] == ", p_formatted))
  }
}, character(1))

p1 <- p1 +
  theme(legend.position = "bottom") + 
  geom_segment(data = sig_bars, 
               aes(x = x, xend = xend, y = y, yend = y),
               size = 0.7) +
  geom_segment(data = sig_bars,
               aes(x = x, xend = x, y = y, yend = y - max_height * 0.02),
               size = 0.7) +
  geom_segment(data = sig_bars,
               aes(x = xend, xend = xend, y = y, yend = y - max_height * 0.02),
               size = 0.7) +
  geom_text(data = sig_bars,
            aes(x = (x + xend) / 2, y = y + max_height * 0.05),
            label = sig_bars$label_text,
            parse = TRUE,
            size = 3.5) +
  coord_cartesian(ylim = c(0, max(sig_bars$y) * 1.05))

print(p1)

ggsave("./Output/Output_figs/K4/K4_SporeSmokeLevel_model_predictions.png", p1, dpi = 600, width = 4.5, height = 4.5)
