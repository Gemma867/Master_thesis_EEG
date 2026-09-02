
# This .R file contains the code of the Q-Q plots from the comparison between models
# with p=7 and fs=1000 Hz (the original parameters).


library(ggplot2)
library(dplyr)

EEGtrial <- readRDS("matrices.rds")

parameters_trials_notfull_p7_phi<-readRDS("parameters_trials_notfull_p7_phi.rds")
parameters_6trials_junt <- readRDS("parameters_6trials_junt.rds")

channels_index <- c(24, 27,  8, 57, 39)
trials <-  c(1, 2, 4, 5, 21, 223)

  

trial_n <- 1:3
p <- 7

qq_data <- list()

for (tn in trial_n) {
  
  eeg_data <- t(EEGtrial[[trials[tn]]])
  
  for (ch in channels_index) {
    
    idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
    
    # Same selected channels for both models
    y <- eeg_data[, idx[1:p]]
    
    # --------------------------------------------------
    # FIRST MODEL
    # --------------------------------------------------
    
    dlmM1 <- buildSignal(
      parameters_6trials_junt[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal <- dropFirst(eegSmo$s)
    residuals <- y[, 1] - signal
    
    shapiro_p <- shapiro.test(residuals)$p.value
    
    # Storage of residuals and Shapiro-Wilk test results for second model
    qq_data[[length(qq_data) + 1]] <- data.frame(
      residuals = residuals,
      channel = ordenSensores$electrode[ch],
      trial = trials[tn],
      model = "First model",
      shapiro_p = shapiro_p,
      shapiro_result = ifelse(
        shapiro_p < 0.05,
        "Does not pass",
        "Passes"
      )
    )
    
    
    # --------------------------------------------------
    # MODIFIED MODEL
    # --------------------------------------------------
    
    dlmM1 <- buildSignalnotfullphineq0(
      parameters_trials_notfull_p7_phi[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal <- dropFirst(eegSmo$s)
    residuals <- y[, 1] - signal
    
    shapiro_p <- shapiro.test(residuals)$p.value
    
    # Storage of residuals and Shapiro-Wilk test results for second model
    qq_data[[length(qq_data) + 1]] <- data.frame(
      residuals = residuals,
      channel = ordenSensores$electrode[ch],
      trial = trials[tn],
      model = "Modified model",
      shapiro_p = shapiro_p,
      shapiro_result = ifelse(
        shapiro_p < 0.05,
        "Does not pass",
        "Passes"
      )
    )
  }
}

qq_df <- bind_rows(qq_data)

qq_df$trial <- factor(
  qq_df$trial,
  levels = trials
)

qq_df$model <- factor(
  qq_df$model,
  levels = c("First model", "Modified model")
)

# Plot
pl2 <- ggplot(
  qq_df,
  aes(
    sample = residuals,
    color = shapiro_result
  )
) +
  
  stat_qq(
    size = 0.7,
    alpha = 0.6
  ) +
  
  stat_qq_line(
    color = "red",
    linewidth = 0.5
  ) +
  
  facet_grid(
    channel ~ trial + model,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "Passes" = "black",
      "Does not pass" = "blue"
    )
  ) +
  
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles",
    color = "Shapiro-Wilk"
  ) +
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw() 

pl2
