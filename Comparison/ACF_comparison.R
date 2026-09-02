
# This .R file contains the code to obtain the ACF plots of the first and second model comparison with p=7 and fs=1000 Hz.

library(ggplot2)
library(dplyr)

EEGtrial <- readRDS("matrices.rds")

parameters_trials_notfull_p7_phi<-readRDS("parameters_trials_notfull_p7_phi.rds")
parameters_6trials_junt <- readRDS("parameters_6trials_junt.rds")

channels_index <- c(24, 27,  8, 57, 39)
trials <-  c(1, 2, 4, 5, 21, 223)

trial_n <- 1:3
p <- 7

acf_data <- list()

for (tn in trial_n) {
  
  eeg_data <- t(EEGtrial[[trials[tn]]])
  
  for (ch in channels_index) {
    
    idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
    
    # Same selected channels for both models
    y <- eeg_data[, idx[1:p]]
    
    # --------------------------------------------------
    # FULL MODEL
    # --------------------------------------------------
    
    dlmM1 <- buildSignal(
      parameters_6trials_junt[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal <- dropFirst(eegSmo$s)
    residuals <- y[, 1] - signal
    
    acf_result <- acf(
      residuals,
      lag.max = 500,
      plot = FALSE
    )
    
    conf <- 2 / sqrt(acf_result$n.used)
    
    acf_data[[length(acf_data) + 1]] <- data.frame(
      lag = acf_result$lag[, 1, 1],
      acf = acf_result$acf[, 1, 1],
      conf_low = -conf,
      conf_high = conf,
      channel = ordenSensores$electrode[ch],
      trial = trials[tn],
      model = "First model"
    )
    
    
    # --------------------------------------------------
    # NOT-FULL MODEL
    # --------------------------------------------------
    
    dlmM1 <- buildSignalnotfullphineq0(
      parameters_trials_notfull_p7_phi[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal <- dropFirst(eegSmo$s)
    residuals <- y[, 1] - signal
    
    acf_result <- acf(
      residuals,
      lag.max = 500,
      plot = FALSE
    )
    
    conf <- 2 / sqrt(acf_result$n.used)
    
    acf_data[[length(acf_data) + 1]] <- data.frame(
      lag = acf_result$lag[, 1, 1],
      acf = acf_result$acf[, 1, 1],
      conf_low = -conf,
      conf_high = conf,
      channel = ordenSensores$electrode[ch],
      trial = trials[tn],
      model = "Modified model"
    )
  }
}

acf_df <- bind_rows(acf_data)

acf_df$trial_model <- paste(
  "Trial", acf_df$trial, acf_df$model
)

trial_model_levels <- as.vector(
  rbind(
    paste("Trial", trials, "First model"),
    paste("Trial", trials, "Modified model")
  )
)

acf_df$trial_model <- factor(
  acf_df$trial_model,
  levels = trial_model_levels
)

acf_df$trial <- factor(acf_df$trial, levels = trials)
acf_df$model <- factor(
  acf_df$model,
  levels = c("First model", "Modified model")
)

pl <- ggplot(acf_df, aes(x = lag, y = acf)) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  
  geom_hline(
    aes(yintercept = conf_high),
    linetype = "dashed",
    color = "blue",
    linewidth = 0.5
  ) +
  
  geom_hline(
    aes(yintercept = conf_low),
    linetype = "dashed",
    color = "blue",
    linewidth = 0.5
  ) +
  
  geom_segment(
    aes(
      xend = lag,
      yend = 0,
      color = model
    ),
    linewidth = 0.3
  ) +
  
  facet_grid(
    channel ~ trial + model,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "First model" = "darkred",
      "Modified model" = "darkgreen"
    )
  ) +
  
  scale_x_continuous(
    breaks = c(0, 250, 500),
    limits = c(0, 500)
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.major = element_line(
      color = "grey70",
      linewidth = 0.4
    ),
    panel.grid.minor = element_line(
      color = "grey85",
      linewidth = 0.25
    ),
    
    axis.text = element_text(size = 9),
    panel.spacing.x = unit(0.8, "lines"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    
    axis.title = element_blank(),
    
    legend.position = "none"
  )
pl
