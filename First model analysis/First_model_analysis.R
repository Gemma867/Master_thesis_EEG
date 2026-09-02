# This .R file contains the code to obtain the plots of the analysis of the Kalman smoother results with the 
# function buildSignal(), p=7 and sampling frequency of 1000 Hz (the original one).

library(GGally)
library(ggplot2)
library(dplyr)
library(dlm)

source("Functions.R")
trials <- c(1,2,4,5,21,223)
EEGtrial<-readRDS("matrices.rds")
parameters_6trials_junt <- readRDS("parameters_6trials_junt.rds")

#### Violin plot of VRR 

# Store results for all trials
results <- vector("list", length(trials))
names(results) <- trials

for (n in seq_along(trials)) {
  
  trial <- trials[n]
  eeg_data <- t(EEGtrial[[trial]])
  
  # Initialize storage for this trial
  var_selected_trial <- vector("list", 64)
  cor_og_trial <- numeric(64)
  varch_trial <- numeric(64)
  var_trial <- numeric(64)
  
  cor_min_trial <- numeric(64)
  var_max_trial <- numeric(64)
  
  cor_mat <- cor(eeg_data)   # compute once per trial
  
  for (ch in 1:64) {
    
    idx <- order(cor_mat[, ch], decreasing = TRUE)
    
    p <- 7
    y <- eeg_data[, idx[1:p]]
    
    dlmM1 <- buildSignal(parameters_6trials_junt[[n]][[ch]])
    eegSmo <- dlmSmooth(y, dlmM1)
    
    senyal <- dropFirst(eegSmo$s)
    
    # Correlation original vs reconstructed
    cor_og_trial[ch] <- cor(senyal, eeg_data[, idx[1]])
    
    # Variance
    var_trial[ch] <- var(eeg_data[, idx[1]])
    
    # Variance change %
    varch_trial[ch] <-
      (var_trial[ch] - var(senyal)) / var_trial[ch] * 100
    
    # Minimum correlation among selected channels
    cor_min_trial[ch] <- min(cor_mat[idx[1:p], ch])
    
    # Variance of selected channels
    var_selected_trial[[ch]] <-
      apply(eeg_data[, idx[1:p]], 2, var)
    
    # Maximum variance excluding target channel
    var_max_trial[ch] <-
      max(apply(eeg_data[, idx[2:p]], 2, var))
  }
  
  # Save everything for this trial
  results[[n]] <- list(
    trial = trial,
    cor_og = cor_og_trial,
    var = var_trial,
    var_change = varch_trial,
    cor_min = cor_min_trial,
    var_selected = var_selected_trial,
    var_max = var_max_trial
  )
}

# Create dataframe
df_all <- bind_rows(
  lapply(results, function(x) {
    data.frame(
      trial = factor(x$trial),
      var_change = x$var_change
    )
  })
)

# Order trials according to input order
df_all$trial <- factor(df_all$trial, levels = trials)

# one violin per trial
ggplot(df_all,
       aes(x = trial,
           y = var_change)) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white"
  ) +
  labs(
    title = "Distribution of variance reduction rate across trials",
    x = "Trial",
    y = "Variance reduction rate (%)"
  ) +
  theme_minimal()


#### Scatter plot VRR against correlation

library(ggplot2)
library(dplyr)
library(ggExtra)

# Collect all values
df_all <- bind_rows(
  lapply(results, function(x) {
    
    keep <- setdiff(1:64, 61)
    
    data.frame(
      vr = x$var_change[keep],
      cor = x$cor_og[keep]
    )
  })
)

# Scatter plot
plot_6<-ggplot(df_all,
               aes(y = cor,
                   x = vr)) +
  
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  
  labs(
    y = "Correlation",
    x = "Variance reduction (%)"
  ) +
  
  theme_minimal()

ggMarginal(plot_6, type="histogram", bins=30,
           yparams = list(fill = "orange"),  # bottom histogram
           xparams = list(fill = "skyblue"))




#### Same plot per trial

library(ggplot2)
library(dplyr)
library(ggrepel)

# Build dataframe
df_cor_vr <- bind_rows(
  lapply(results, function(x) {
    
    keep <- setdiff(1:64, 61)
    
    df <- data.frame(
      trial = factor(x$trial),
      electrode = keep,
      cor = x$cor_og[keep],
      vr = x$var_change[keep]
    )
    
    # Labels only for selected points
    df$label <- ifelse(
      df$cor < 0.75 | (df$vr < 10 | df$vr > 60),
      ordenSensores$electrode,
      NA
    )
    
    df
  })
)

# One figure with one panel per trial
ggplot(df_cor_vr,
       aes(x = vr,
           y = cor, color=trial)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    aes(label = label),
    na.rm = TRUE,
    size = 3,
    show.legend = FALSE
  ) +
  
  facet_wrap(~ trial) +
  labs(
    y = "Correlation",
    x = "Variance reduction rate (%)"
  ) +
  theme_minimal()+
  theme(
    legend.position = "none"
  )


#### GGpairs plot

library(GGally)
library(dplyr)

# Build dataframe
df_pairs <- bind_rows(
  lapply(results, function(x) {
    
    keep <- setdiff(1:64, 61)
    
    data.frame(
      VR = x$var_change[keep],
      cor = x$cor_og[keep],
      LowestCor = x$cor_min[keep],
      Variance = x$var[keep]
    )
  })
)

# Pair plot
ggpairs(
  df_pairs,
  columns = c("VR", "cor", "LowestCor", "Variance"),
  upper = list(
    continuous = wrap(
      "cor",
      size = 4
    )
  ),
  lower = list(
    continuous = wrap(
      "points",
      alpha = 0.5,
      size = 1
    )
  ),
  diag = list(
    continuous = wrap(
      "densityDiag",
      alpha = 0.6
    )
  )
)


##### Barplot of correlations between residuals of signal from the same model

library(dplyr)
library(tidyr)
library(ggplot2)

# Initialize storage of innovations and correlation between reference innovations and other surrounding signals from the same model
innov_trial <- vector("list", 64)
corr_innov_trial <- vector("list", 64)
results_res<-vector("list",6)
trials<-c(1,2,4,5,21,223)
for (n in seq_along(trials)){
  for (ch in 1:64) {
    trial <- trials[n]
    eeg_data <- t(EEGtrial[[trial]])
    idx <- order(cor(eeg_data)[,ch], decreasing = TRUE)
    
    p <- 7
    y <- eeg_data[, idx[1:p]]
    
    dlmM1 <- buildSignal(parameters_6trials_junt[[n]][[ch]])
    
    # Filter
    eegFilt <- dlmSmooth(y, dlmM1)
    
    innov <- y-unlist(eegFilt)
    innov <- data.matrix(innov)
    # Store innovations
    innov_trial[[ch]] <- innov
    
    # Correlation between innovations of channels 
    corr_innov_trial[[ch]] <-
      cor(innov, use = "pairwise.complete.obs")
  }
  
  
  # Save
  results_res[[n]] <- list(
    trial = trial,
    innovations = innov_trial,
    corr_innov = corr_innov_trial
  )}

df<-data.frame(Trial=rep(1:6,each=64*6),Ref=rep(rep(ordenSensores$electrode,each=6),times=6),Top6=rep(1,times=6*64*6),Cor_Top6=rep(1,times=6*64*6))

for (n in 1:6){
  trial<-trials[n]
  eeg_data <- t(EEGtrial[[trial]])
  for (i in 1:64){
    for (j in 1:6){
      idx <- order(cor(eeg_data)[,i], decreasing = TRUE)
      df[(i-1)*6+j+(n-1)*64*6,4]<-results_res[[n]]$corr_innov[[i]][1+j,1]
      df[(i-1)*6+j+(n-1)*64*6,3]<-ordenSensores$electrode[idx[1+j]]
    }
  }
  print(trial)
}

# Count occurrences above thresholds
plot_df <- df %>%
  filter(!is.na(Ref)) %>% # Not counting the reference channel
  group_by(Ref) %>%
  summarise(
    `0.3-0.4` = sum(Cor_Top6 >= 0.3 & Cor_Top6 < 0.4),
    `0.4-0.5` = sum(Cor_Top6 >= 0.4 & Cor_Top6 < 0.5),
    `0.5-0.6` = sum(Cor_Top6 >= 0.5 & Cor_Top6 < 0.6),
    `0.6<` = sum(Cor_Top6 >= 0.6)
  ) %>%
  pivot_longer(
    cols = starts_with("0"),
    names_to = "Threshold",
    values_to = "Count"
  )

# Plot
ggplot(plot_df, aes(x = Ref, y = Count, fill = Threshold)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(
      "0.3-0.4" = "skyblue",
      "0.4-0.5" = "orange",
      "0.5-0.6" = "red",
      "0.6<" = "darkred"
    )
  ) +
  labs(
    x = "Reference Channel",
    y = "Number of correlation values in the ranges",
    fill = "Correlation of residuals \nof top 6 channels with reference \nchannel"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) 

#### 2D layout maximum correlation achieved between residuals of the reference channel and residuals of other (p-1) channels from the model

library(dplyr)

df <- df %>%
  left_join(chan_coords, by = c("Ref" = "electrode"))

df_channels <- df %>%
  group_by(Ref, x, y) %>%
  summarise(
    max_cor = max(Cor_Top6),
    .groups = "drop"
  ) %>%
  mutate(
    Group = case_when(
      max_cor > 0.6 ~ ">0.6",
      max_cor > 0.5 ~ "0.5-0.6",
      max_cor > 0.4 ~ "0.4–0.5",
      max_cor > 0.3 ~ "0.3–0.4",
      TRUE ~ "<0.3"
    )
  )

library(ggplot2)

ggplot(df_channels, aes(x = x, y = y, color = Group)) +
  geom_point(size = 4) +
  coord_equal() +
  scale_color_manual(values = c(
    ">0.6" = "darkred",
    "0.5-0.6" = "red",
    "0.4–0.5" = "orange",
    "0.3–0.4" = "skyblue",
    "<0.3" = "grey80"
  )) +
  theme_void() +
  geom_text_repel(
    aes(label = Ref),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.6,
    max.overlaps = Inf,
    color = "black"
  ) +
  labs(
    color = "Max correlation"
  )+
  theme(
    legend.text = element_text(size = 10)
  )


#### ACF

chs <- seq(1, 48, by = 3)   # 1, 4, 7, ..., 58 (20 channels)

par(mfrow = c(4, 4),
    mar = c(3, 3, 3, 1),  
    cex.main = 1.6)         

trial <- 1 # trial to plot
eeg_data <- t(EEGtrial[[trial]])
p <- 7 # p is always 7 in this first model

for (ch in chs) {
  
  idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
  y <- eeg_data[, idx[1:p]]
  
  dlmM1 <- buildSignal(parameters_6trials_junt[[trial]][[ch]])
  eegSmo <- dlmSmooth(y, dlmM1)
  senyal <- dropFirst(eegSmo$s)
  residuals <- y[,1]-senyal
  
  acf(residuals,
      lag.max = 500,
      main = paste0(as.character(ordenSensores$electrode[ch])),
      cex.main = 1.6)
}

#### PSD


library(ggplot2)
library(patchwork)

chs <- seq(25, 48, by = 3)

trial <- 1
eeg_data <- t(matrices[[trial]])
p <- 7

Fs <- 1000  # Sampling frequency

plots <- vector("list", length(chs))

for (i in seq_along(chs)) {
  
  ch <- chs[i]
  
  idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
  y <- eeg_data[, idx[1:p]]
  
  dlmM1 <- buildSignal(parameters_6trials_junt[[trial]][[ch]])
  eegSmo <- dlmSmooth(y, dlmM1)
  
  original <- y[, 1]
  denoised <- dropFirst(eegSmo$s)
  
  # Spectrum of original signal
  sp_orig <- spectrum(original, log="no", spans = c(3, 3), plot = FALSE, detrend = TRUE)
  # Spectrum of denoised signal
  sp_den  <- spectrum(denoised, log="no", spans = c(3, 3), plot = FALSE, detrend = TRUE)
  
  # Store psd and corresponding frequency in dataframes
  df_orig <- data.frame(
    Frequency = sp_orig$freq * Fs,
    PSD = sp_orig$spec
  )
  
  df_den <- data.frame(
    Frequency = sp_den$freq * Fs,
    PSD = sp_den$spec
  )
  
  # plots is a list of plots where each plot contains the spectrum of the original and denoised superimposed
  plots[[i]] <- ggplot() +
    geom_line(data = df_orig,
              aes(Frequency, PSD),
              colour = "black",
              linewidth = 0.8) +
    geom_line(data = df_den,
              aes(Frequency, PSD),
              colour = "red",
              linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 60),breaks = seq(0, 300, by = 10)) +
    labs(
      title = ordenSensores$electrode[ch],
      x = "Frequency (Hz)",
      y = "PSD"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        vjust = 1
      )
    )
}

wrap_plots(plots, ncol = 4)

#### Histogram power reduction

# Definition of upper and lower limits bands
Fs=1000 # sampling frequency
bands <- data.frame(
  band = c("Delta", "Theta", "Alpha", "Beta"),
  low  = c(0.5, 4, 8, 13)/Fs,
  high = c(4,   8, 13, 30)/Fs
)

# Function to compute the power in one frequency band
band_power <- function(sp, f_low, f_high) {
  idx <- sp$freq >= f_low & sp$freq < f_high
  sum(sp$spec[idx])
}

# Function to compute the total power of the spectrum
total_power <- function(sp, f_low, f_high) {
  sum(sp$spec)
}

chs <- 1:64
trials<-c(1,2,4,5,21,223)

results_freq_bands <- vector("list", length(trials))

for (n in seq_along(trials)) {
  
  trial <- trials[n]
  eeg_data <- t(matrices[[trial]])
  trial_results <- data.frame()
  
  for (ch in 1:64) {
    
    idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
    y <- eeg_data[, idx[1:p]]
    
    dlmM1 <- buildSignal(parameters_6trials_junt[[n]][[ch]])
    eegSmo <- dlmSmooth(y, dlmM1)
    
    original <- y[,1]
    denoised <- dropFirst(eegSmo$s)
    
    # Spectrum of original signal
    sp_orig <- spectrum(original, spans = c(3,3),
                        plot = FALSE, detrend = TRUE)
    
    # Spectrum of denoised signal
    sp_den <- spectrum(denoised, spans = c(3,3),
                       plot = FALSE, detrend = TRUE)
    
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
    den_total  <- total_power(sp_den)
    
    # Storage of power reduction (change_pct) and other quantities
    trial_results <- rbind(
      trial_results,
      data.frame(
        trial = trial,
        channel = ch,
        band = bands$band,
        original = orig,
        denoised = den,
        rel_orig = orig / orig_total,
        rel_denoised = den / den_total,
        change_pct = 100 * (orig-den) / orig # Power reduction
      )
    )
  }
  
  results_freq_bands[[n]] <- trial_results
}

## Per waveband

ggplot(df_plot, aes(change_pct)) +
  geom_histogram(bins = 30, color = "black", fill = "darkolivegreen3") +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_x_continuous(breaks = seq(-100, 100, by = 20))+
  facet_grid(cols = vars(band))+
  theme_minimal() +
  labs(
    x = "Change in band power (%)",
    y = "Count"
  )


## Per trial and waveband

library(dplyr)
library(purrr)
library(ggplot2)

df_plot <- bind_rows(results_freq_bands)

df_plot$band <- factor(
  df_plot$band,
  levels = c("Delta", "Theta", "Alpha", "Beta")
)

df_plot <- df_plot %>% filter(band %in% c("Delta", "Theta", "Alpha", "Beta"))

medians <- df_plot %>%
  group_by(trial, band) %>%
  summarise(
    median = median(change_pct, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(df_plot, aes(change_pct)) +
  geom_histogram(bins = 30, color = "black", fill = "coral") +
  geom_vline(
    data = medians,
    aes(xintercept = median),
    color = "red",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_x_continuous(breaks = seq(-100, 100, by = 20))+
  facet_grid(rows = vars(trial), cols = vars(band))+
  theme_minimal() +
  labs(
    x = "Change in band power (%)",
    y = "Count"
  )


#### Plot of signals in black and denoised signal in red

par(mfrow = c(6, 2), mar=c(2,2,1,1), oma = c(0, 0, 5, 1)) # outer margin for global title
for (i in idx[1:p]){
  plot(ts(eeg_data[,i]),
       main=paste0(i,"-",ordenSensores$electrode[i]),
       ylim=c(-20,20))
  lines(senyal,col=2)
}
mtext(" Denoised target Signal (red) and most correlated EEG channels",
      outer = TRUE, line=2, cex = 1.4, font = 1)
