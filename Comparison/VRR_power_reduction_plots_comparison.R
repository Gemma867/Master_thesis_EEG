
library(ggplot2)
library(dplyr)
library(dlm)

EEGtrial <- readRDS("matrices.rds")
source("Functions.R")

parameters_trials_notfull_p7_phi<-readRDS("parameters/parameters_trials_notfull_p7_phi.rds")
parameters_6trials_junt <- readRDS("parameters/parameters_6trials_junt.rds")

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

### Boxplot band power reduction comparison between first and second model

p <- 7
Fs <- 1000 # sampling frequency

# Models and corresponding build functions
models <- list(
  "First model" = list(
    pars = parameters_6trials_junt,
    build = buildSignal
  ),
  "Modified model" = list(
    pars = parameters_trials_notfull_p7_phi,
    build = buildSignalnotfullphineq0
  )
)

# Frequency bands
bands <- data.frame(
  band = c("Delta", "Theta", "Alpha", "Beta"),
  low  = c(0.5, 4, 8, 13) / Fs,
  high = c(4, 8, 13, 30) / Fs
)

results_freq_bands <- list()

for (model_name in names(models)) {
  
  model <- models[[model_name]]
  
  for (n in seq_along(trials)) {
    
    trial <- trials[n]
    eeg_data <- t(EEGtrial[[trial]])
    
    for (ch in channels_index) {
      
      idx <- order(
        cor(eeg_data)[, ch],
        decreasing = TRUE
      )
      
      # Select p channels
      y <- eeg_data[, idx[1:p], drop = FALSE]
      
      # Reference channel only for p = 1
      if (p == 1) {
        y <- matrix(y[, 1], ncol = 1)
      }
      
      # Build model
      dlmM1 <- model$build(
        model$pars[[n]][[ch]]
      )
      
      # Smooth
      eegSmo <- dlmSmooth(y, dlmM1)
      
      # Reference channel
      if (p == 1) {
        original <- y[, 1]
      } else {
        original <- y[, 1]
      }
      
      # Denoised signal
      denoised <- dropFirst(eegSmo$s)
      
      # Spectrum
      sp_orig <- spectrum(
        original,
        plot = FALSE,
        detrend = TRUE
      )
      
      sp_den <- spectrum(
        denoised,
        plot = FALSE,
        detrend = TRUE
      )
      
      # Band powers
      orig <- mapply(
        band_power,
        MoreArgs = list(sp = sp_orig),
        bands$low,
        bands$high
      )
      
      den <- mapply(
        band_power,
        MoreArgs = list(sp = sp_den),
        bands$low,
        bands$high
      )
      
      orig_total <- total_power(sp_orig)
      den_total <- total_power(sp_den)
      
      # Store
      results_freq_bands[[length(results_freq_bands) + 1]] <-
        data.frame(
          trial = trial,
          channel = ch,
          Fs = Fs,
          model = model_name,
          band = bands$band,
          original = orig,
          denoised = den,
          rel_orig = orig / orig_total,
          rel_denoised = den / den_total,
          change_pct = 100 * (orig - den) / orig
        )
    }
  }
}

df_plot <- bind_rows(results_freq_bands)

df_plot$trial <- factor(
  df_plot$trial,
  levels = c(1, 2, 4, 5, 21, 223)
)

df_plot$band <- factor(
  df_plot$band,
  levels = c("Delta", "Theta", "Alpha", "Beta"),
  labels = c(
    "Delta (0.5–4 Hz)",
    "Theta (4–8 Hz)",
    "Alpha (8–13 Hz)",
    "Beta (13–30 Hz)"
  )
)

df_plot$model <- factor(
  df_plot$model,
  levels = c("First model", "Modified model")
)

# Plot
ggplot(
  df_plot,
  aes(
    x = model,
    y = change_pct,
    fill = model
  )
) +
  geom_boxplot(
    #alpha = 0.6,
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  facet_grid(
    cols = vars(band)
  ) +
  scale_y_continuous(
    limits = c(-20, 75),
    breaks = seq(-20, 75, by = 10)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    legend.position = "none"
  ) +
  labs(
    x = "",
    y = "Band power reduction (%)"
  )

 
