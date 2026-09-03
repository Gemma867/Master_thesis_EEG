
library(ggplot2)
library(dplyr)
library(dlm)

EEGtrial <- readRDS("matrices.rds")

parameters_trials_notfull_p7_phi<-readRDS("parameters_trials_notfull_p7_phi.rds")
parameters_6trials_junt <- readRDS("parameters_6trials_junt.rds")

channels_index <- c(24, 27,  8, 57, 39)
trials <-  c(1, 2, 4, 5, 21, 223)

trial_n <- 1:6
p <- 7

vrr_data <- list()

for (tn in trial_n) {
  
  eeg_data <- t(EEGtrial[[trials[tn]]])
  
  for (ch in channels_index) {
    
    idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
    
    y <- eeg_data[, idx[1:p]]
    
    # ==================================================
    # FULL MODEL
    # ==================================================
    
    dlmM1 <- buildSignal(
      parameters_6trials_junt[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal_full <- dropFirst(eegSmo$s)
    
    vrr_full <- (
      var(y[, 1]) - var(signal_full)
    ) / var(y[, 1]) * 100
    
    
    # ==================================================
    # NOT-FULL MODEL
    # ==================================================
    
    dlmM1 <- buildSignalnotfullphineq0(
      parameters_trials_notfull_p7_phi[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal_notfull <- dropFirst(eegSmo$s)
    
    vrr_notfull <- (
      var(y[, 1]) - var(signal_notfull)
    ) / var(y[, 1]) * 100
    
    
    # Store
    vrr_data[[length(vrr_data) + 1]] <- data.frame(
      trial = trials[tn],
      channel = ordenSensores$electrode[ch],
      VRR_full = vrr_full,
      VRR_notfull = vrr_notfull
    )
  }
}

vrr_df <- bind_rows(vrr_data)

#### Scatter VRR plot

p1<-ggplot(vrr_df, aes(x = VRR_full, y = VRR_notfull)) +
  
  geom_point(
    size = 2,
    alpha = 0.6,
    color = "blue"
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  labs(
    x = "VRR — First model (%)",
    y = "VRR — Modified model (%)"
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
    )
  )

p1

#### Scatter correlation plot

cor_data <- list()

for (tn in 1:6) {
  
  eeg_data <- t(EEGtrial[[trials[tn]]])
  
  for (ch in channels_index) {
    
    idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
    
    y <- eeg_data[, idx[1:p]]
    
    # ==================================================
    # FULL MODEL
    # ==================================================
    
    dlmM1 <- buildSignal(
      parameters_6trials_junt[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal_full <- dropFirst(eegSmo$s)
    
    cor_full <- cor(
      y[, 1],
      signal_full
    )
    
    
    # ==================================================
    # NOT-FULL MODEL
    # ==================================================
    
    dlmM1 <- buildSignalnotfullphineq0(
      parameters_trials_notfull_p7_phi[[tn]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    
    signal_notfull <- dropFirst(eegSmo$s)
    
    cor_notfull <- cor(
      y[, 1],
      signal_notfull
    )
    
    
    # Store
    cor_data[[length(cor_data) + 1]] <- data.frame(
      trial = trials[tn],
      channel = ordenSensores$electrode[ch],
      Cor_full = cor_full,
      Cor_notfull = cor_notfull
    )
  }
}

cor_df <- bind_rows(cor_data)

p2<-ggplot(
  cor_df,
  aes(x = Cor_full, y = Cor_notfull)
) +
  
  geom_point(
    size = 2,
    alpha = 0.6,
    color = "blue"
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  labs(
    x = "Correlation — First model",
    y = "Correlation — Modified model"
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
    )
  )

p2
 