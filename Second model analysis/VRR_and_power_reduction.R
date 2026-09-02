
# This .R shows the code to compute the VRR and power reduction metrics. 
# For the VRR, one figure with boxplots in function of the sampling frequency (fs) is coded. The p and list of models must be changed accordingly.
# For the power reduction, two figures (histograms and boxplots in function of the wavebands) are coded. The Fs, p and list of models must be changed accordingly.

# Loading of parameters estimations.

library(dplyr)
library(ggplot2)
library(gsignal)

EEGtrial <- readRDS("matrices.rds")

parameters_trials_notfull_p3_phi<-readRDS("parameters_trials_notfull_p3_phi.rds")
parameters_trials_notfull_p5_phi<-readRDS("parameters_trials_notfull_p5_phi.rds")
parameters_trials_notfull_p7_phi<-readRDS("parameters_trials_notfull_p7_phi.rds")
parameters_trials_notfull_p1_phi<-readRDS("parameters_trials_notfull_p1_phi.rds")

parameters_trials_notfull_p3_phi_500<-readRDS("parameters_trials_notfull_p3_phi_500.rds")
parameters_trials_notfull_p5_phi_500<-readRDS("parameters_trials_notfull_p5_phi_500.rds")
parameters_trials_notfull_p7_phi_500<-readRDS("parameters_trials_notfull_p7_phi_500.rds")
parameters_trials_notfull_p1_phi_500<-readRDS("parameters_trials_notfull_p1_phi_500.rds")

parameters_trials_notfull_p3_phi_250<-readRDS("parameters_trials_notfull_p3_phi_250.rds")
parameters_trials_notfull_p5_phi_250<-readRDS("parameters_trials_notfull_p5_phi_250.rds")
parameters_trials_notfull_p7_phi_250<-readRDS("parameters_trials_notfull_p7_phi_250.rds")
parameters_trials_notfull_p1_phi_250<-readRDS("parameters_trials_notfull_p1_phi_250.rds")

channels_index <- c(24, 27,  8, 57, 39)
trials <-  c(1, 2, 4, 5, 21, 223)

#Computation of VRR

#Example for p = 3

p <- 3

if (p == 1) {
  
  models <- list(
    list(
      ds = 1,
      pars = parameters_trials_notfull_p1_phi,
      title = "1000 Hz"
    ),
    list(
      ds = 2,
      pars = parameters_trials_notfull_p1_phi_500,
      title = "500 Hz"
    ),
    list(
      ds = 4,
      pars = parameters_trials_notfull_p1_phi_250,
      title = "250 Hz"
    )
  )
  
} else if (p == 3) {
  
  models <- list(
    list(
      ds = 1,
      pars = parameters_trials_notfull_p3_phi,
      title = "1000 Hz"
    ),
    list(
      ds = 2,
      pars = parameters_trials_notfull_p3_phi_500,
      title = "500 Hz"
    ),
    list(
      ds = 4,
      pars = parameters_trials_notfull_p3_phi_250,
      title = "250 Hz"
    )
  )
  
} else if (p == 5) {
  
  models <- list(
    list(
      ds = 1,
      pars = parameters_trials_notfull_p5_phi,
      title = "1000 Hz"
    ),
    list(
      ds = 2,
      pars = parameters_trials_notfull_p5_phi_500,
      title = "500 Hz"
    ),
    list(
      ds = 4,
      pars = parameters_trials_notfull_p5_phi_250,
      title = "250 Hz"
    )
  )
  
} else if (p == 7) {
  
  models <- list(
    list(
      ds = 1,
      pars = parameters_trials_notfull_p7_phi,
      title = "1000 Hz"
    ),
    list(
      ds = 2,
      pars = parameters_trials_notfull_p7_phi_500,
      title = "500 Hz"
    ),
    list(
      ds = 4,
      pars = parameters_trials_notfull_p7_phi_250,
      title = "250 Hz"
    )
  )
}

vrr_results <- list()

for (trial_n in seq_along(trials)) {
  
  eeg_data <- t(EEGtrial[[trials[trial_n]]])
  cor_mat <- cor(eeg_data)
  
  for (ch in channels_index) {
    
    idx <- order(
      cor_mat[, ch],
      decreasing = TRUE
    )
    
    for (m in models) {
      
      y <- eeg_data[, idx[1:p]]
      
      # Downsample
      if (p==1){
      y <- downsample(y, m$ds)}
      if (p>1){
        y <- apply(y, 2, downsample, m$ds)
      }
      
      # Model parameters
      par <- m$pars[[trial_n]][[ch]]
      
      # Build model
      if (p==1){
      dlmM1 <- buildSignalnotfullphineq0_p1(par)}
      if (p>1){
      dlmM1 <- buildSignalnotfullphineq0(par)}
      
      # Smooth signal
      eegSmo <- dlmSmooth(y, dlmM1)
      
      # Estimated signal
      signal <- dropFirst(eegSmo$s)
      
      # Residuals
      if (p==1){
      residuals <- y - signal
      var_og <- var(y)}
      if (p>1){
      residuals <- y[,1] - signal
      var_og <- var(y[,1])}
      
      # Variances
      var_res <- var(residuals)
      var_den <- var(signal)
      
      # Store results
      vrr_results[[length(vrr_results) + 1]] <- data.frame(
        trial = trials[trial_n],
        channel = ordenSensores$electrode[ch],
        fs = m$title,
        p = 1,
        var_res = var_res, # variance of residuals
        var_og = var_og, # variance of original reference signal
        var_den = var_den, # variance of denoised signal
        vRR_res = var_res / var_og * 100, # fraction of original variance represented by residual variance 
        VRR = (var_og - var_den) / var_og * 100 # variance reduction
      )
    }
  }
}

vrr_df <- bind_rows(vrr_results)

# Plot of the results

vrr_df$fs <- factor(vrr_df$fs, levels = c("1000 Hz", "500 Hz", "250 Hz"))
ggplot(vrr_df, aes(x = fs, y = VRR, fill = fs)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "Variance Reduction Rate (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )+
  labs(
    title = paste0("p=",p)
  )


# Frequency bands power reduction

p <- 3

Fs_values <- c(1000, 500, 250)

band_power <- function(sp, f_low, f_high) {
  idx <- sp$freq >= f_low & sp$freq < f_high
  sum(sp$spec[idx])
}

total_power <- function(sp) {
  sum(sp$spec)
}

results_freq_bands <- list()

for (Fs in Fs_values) {
  
  # Downsampling factor
  ds <- 1000 / Fs
  
  # Frequency bands
  bands <- data.frame(
    band = c("Delta", "Theta", "Alpha", "Beta"),
    low  = c(0.5, 4, 8, 13) / Fs,
    high = c(4, 8, 13, 30) / Fs
  )
  
  # Select parameters according to p and Fs
  if (p == 1) {
    pars <- switch(
      as.character(Fs),
      "1000" = parameters_trials_notfull_p1_phi,
      "500"  = parameters_trials_notfull_p1_phi_500,
      "250"  = parameters_trials_notfull_p1_phi_250
    )
  } 
  if (p==3){
    pars <- switch(
      as.character(Fs),
      "1000" = parameters_trials_notfull_p3_phi,
      "500"  = parameters_trials_notfull_p3_phi_500,
      "250"  = parameters_trials_notfull_p3_phi_250
    )
  }
  if (p==5){
    pars <- switch(
      as.character(Fs),
      "1000" = parameters_trials_notfull_p5_phi,
      "500"  = parameters_trials_notfull_p5_phi_500,
      "250"  = parameters_trials_notfull_p5_phi_250
    )
  }
  if (p==7){
    pars <- switch(
      as.character(Fs),
      "1000" = parameters_trials_notfull_p7_phi,
      "500"  = parameters_trials_notfull_p7_phi_500,
      "250"  = parameters_trials_notfull_p7_phi_250
    )
  }
  
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
      
      # Downsample
      if (p==1){
        y<-downsample(y,ds)
      } else {
      y <- apply(
        y,
        2,
        downsample,
        ds
      )}
      
      # Build model
      if (p == 1) {
        dlmM1 <- buildSignalnotfullphineq0_p1(
          pars[[n]][[ch]]
        )
      } else {
        dlmM1 <- buildSignalnotfullphineq0(
          pars[[n]][[ch]]
        )
      }
      
      # Smooth
      eegSmo <- dlmSmooth(y, dlmM1)
      
      # Reference channel
      if (p==1){
        original<-y
      } else {
      original <- y[, 1]}
      
      # Denoised signal
      denoised <- dropFirst(eegSmo$s)
      
      # Spectrum
      sp_orig <- spectrum(
        original,
        spans = c(3, 3),
        plot = FALSE,
        detrend = TRUE
      )
      
      sp_den <- spectrum(
        denoised,
        spans = c(3, 3),
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

# Histogram

ggplot(df_plot, aes(change_pct)) +
  geom_histogram(bins = 30, color = "black", fill = "darkolivegreen3") +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  #scale_x_continuous(breaks = seq(-100, 100, by = 20))+
  facet_grid(cols = vars(band))+
  theme_minimal() +
  labs(
    x = "Change in band power (%)",
    y = "Count",
    title = paste0("p=",p)
  )

# Boxplot

df_plot$Fs <- factor(
  df_plot$Fs,
  levels = c(1000, 500, 250)
)

ggplot(
  df_plot,
  aes(
    x = Fs,
    y = change_pct,
    fill = Fs
  )
) +
  geom_boxplot(
    alpha = 0.6,
    color = "black"
  ) +
  facet_grid(
    cols = vars(band)
  ) +
  scale_y_continuous(
    limits = c(-50, 75),
    breaks = seq(-50, 75, by = 25)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    legend.position = "none"
  ) +
  labs(
    x = "",
    y = "Change in band power (%)",
    title = paste0("p = ", p)
  )

