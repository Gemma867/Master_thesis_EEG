
# This .R file contains the code to compute the ACF function once the parameters are estimated, as well as the Ljung-Box test.

# Loading of parameters estimations.


library(ggplot2)
library(dplyr)

source("Functions.R")

EEGtrial <- readRDS("matrices.rds")

parameters_trials_notfull_p3_phi<-readRDS("parameters/parameters_trials_notfull_p3_phi.rds")
parameters_trials_notfull_p5_phi<-readRDS("parameters/parameters_trials_notfull_p5_phi.rds")
parameters_trials_notfull_p7_phi<-readRDS("parameters/parameters_trials_notfull_p7_phi.rds")
parameters_trials_notfull_p1_phi<-readRDS("parameters/parameters_trials_notfull_p1_phi.rds")

parameters_trials_notfull_p3_phi_500<-readRDS("parameters/parameters_trials_notfull_p3_phi_500.rds")
parameters_trials_notfull_p5_phi_500<-readRDS("parameters/parameters_trials_notfull_p5_phi_500.rds")
parameters_trials_notfull_p7_phi_500<-readRDS("parameters/parameters_trials_notfull_p7_phi_500.rds")
parameters_trials_notfull_p1_phi_500<-readRDS("parameters/parameters_trials_notfull_p1_phi_500.rds")

parameters_trials_notfull_p3_phi_250<-readRDS("parameters/parameters_trials_notfull_p3_phi_250.rds")
parameters_trials_notfull_p5_phi_250<-readRDS("parameters/parameters_trials_notfull_p5_phi_250.rds")
parameters_trials_notfull_p7_phi_250<-readRDS("parameters/parameters_trials_notfull_p7_phi_250.rds")
parameters_trials_notfull_p1_phi_250<-readRDS("parameters/parameters_trials_notfull_p1_phi_250.rds")

channels_index <- c(24, 27,  8, 57, 39)
trials <-  c(1, 2, 4, 5, 23, 221)

# ACF subplots for trial_n (to choose from 1 to 6) and p = 1 for all sampling frequencies.


trial_n <- 6
p <- 1

models <- list(
  list(ds = 1,
       pars = parameters_trials_notfull_p1_phi,
       title = "1000 Hz"),
  list(ds = 2,
       pars = parameters_trials_notfull_p1_phi_500,
       title = "500 Hz"),
  list(ds = 4,
       pars = parameters_trials_notfull_p1_phi_250,
       title = "250 Hz")
)

eeg_data <- t(EEGtrial[[trials[trial_n]]]) # EEG data from trial selected

acf_data <- list()

for (ch in channels_index) { # loop across the 5 channels
  
  idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
  
  for (m in models) { # loop across the 3 models (different fs)
    
    y <- eeg_data[, idx[1:p]]
    y <- downsample(y, m$ds) # Downsample of the signal by the factor ds of the model
    
    # Build the filtered signal with the corresponding parameters
    dlmM1 <- buildSignalnotfullphineq0_p1(
      m$pars[[trial_n]][[ch]]
    )
    
    # Smooth the signal
    eegSmo <- dlmSmooth(y, dlmM1)
    
    # Obtain the reference smoothed signal
    signal <- dropFirst(eegSmo$s)
    
    # Compute residuals
    residuals <- y - signal
    
    acf_result <- acf(
      residuals,
      lag.max = 500,
      plot = FALSE
    )
    
    # Extract confidence interval from acf object
    conf <- 1.96 / sqrt(acf_result$n.used)
    
    # Store acf results in a dataframe
    acf_data[[length(acf_data) + 1]] <- data.frame(
      lag = acf_result$lag[, 1, 1],
      acf = acf_result$acf[, 1, 1],
      conf_low = -conf,
      conf_high = conf,
      channel = ordenSensores$electrode[ch],
      ds = m$title
    )
  }
}

acf_df <- bind_rows(acf_data)

acf_df$channel <- factor(
  acf_df$channel,
  levels = ordenSensores$electrode[channels_index]
)

acf_df$ds <- factor(
  acf_df$ds,
  levels = c("1000 Hz", "500 Hz", "250 Hz")
)

# Plot of the acf
pl<-ggplot(acf_df, aes(x = lag, y = acf)) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  
  geom_hline(
    aes(yintercept = conf_high),
    linetype = "dashed",
    color="blue",
    linewidth = 0.5
  ) +
  
  geom_hline(
    aes(yintercept = conf_low),
    linetype = "dashed",
    color="blue",
    linewidth = 0.5
  ) +
  
  geom_segment(
    aes(xend = lag, yend = 0),
    linewidth = 0.3
  ) +
  
  facet_grid(
    channel ~ ds,
    scales = "free_y"
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
    panel.spacing = unit(0.2, "lines"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.title = element_blank()
  )+
  labs(
    title = "Autocorrelation Function of residuals with p=1 and trial=223"
  ) 

pl

# ACF subplots for trial_n (to choose from 1 to 6) and p > 1 for all sampling frequencies. The code is very similar to the p = 1 case, except for changes in the downsampling and computation of residuals.

trial_n <- 2
p <- 5

# for p = 5, 7, the parameters list must be changed accordingly
models <- list(
  list(ds = 1,
       pars = parameters_trials_notfull_p5_phi,
       title = "1000 Hz"),
  list(ds = 2,
       pars = parameters_trials_notfull_p5_phi_500,
       title = "500 Hz"),
  list(ds = 4,
       pars = parameters_trials_notfull_p5_phi_250,
       title = "250 Hz")
)

eeg_data <- t(EEGtrial[[trials[trial_n]]])

acf_data <- list()

for (ch in channels_index) {
  
  idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
  
  for (m in models) {
    
    y <- eeg_data[, idx[1:p]]
    y <- apply(y, 2, function(x) downsample(x, m$ds)) # change with respect p = 1
    
    dlmM1 <- buildSignalnotfullphineq0(
      m$pars[[trial_n]][[ch]]
    )
    
    eegSmo <- dlmSmooth(y, dlmM1)
    signal <- dropFirst(eegSmo$s)
    residuals <- y[, 1] - signal # change with respect p = 1
    
    acf_result <- acf(
      residuals,
      lag.max = 500,
      plot = FALSE
    )
    
    conf <- 1.96 / sqrt(acf_result$n.used)
    
    acf_data[[length(acf_data) + 1]] <- data.frame(
      lag = acf_result$lag[, 1, 1],
      acf = acf_result$acf[, 1, 1],
      conf_low = -conf,
      conf_high = conf,
      channel = ordenSensores$electrode[ch],
      ds = m$title
    )
  }
}

acf_df <- bind_rows(acf_data)

acf_df$channel <- factor(
  acf_df$channel,
  levels = ordenSensores$electrode[channels_index]
)

acf_df$ds <- factor(
  acf_df$ds,
  levels = c("1000 Hz", "500 Hz", "250 Hz")
)

pl<-ggplot(acf_df, aes(x = lag, y = acf)) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  
  geom_hline(
    aes(yintercept = conf_high),
    linetype = "dashed",
    color="blue",
    linewidth = 0.5
  ) +
  
  geom_hline(
    aes(yintercept = conf_low),
    linetype = "dashed",
    color="blue",
    linewidth = 0.5
  ) +
  
  geom_segment(
    aes(xend = lag, yend = 0),
    linewidth = 0.3
  ) +
  
  facet_grid(
    channel ~ ds,
    scales = "free_y"
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
    panel.spacing = unit(0.2, "lines"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.title = element_blank()
  )+
  labs(
    title = "Autocorrelation Function of residuals with p=5 and trial=2"
  ) 

pl

#### Ljung-Box test

p_values <- c(1, 3, 5, 7)

# All parameters are grouped in a list.
models <- list(
  list(
    ds = 1,
    pars = list(
      parameters_trials_notfull_p1_phi,
      parameters_trials_notfull_p3_phi,
      parameters_trials_notfull_p5_phi,
      parameters_trials_notfull_p7_phi
    ),
    title = "1000 Hz"
  ),
  list(
    ds = 2,
    pars = list(
      parameters_trials_notfull_p1_phi_500,
      parameters_trials_notfull_p3_phi_500,
      parameters_trials_notfull_p5_phi_500,
      parameters_trials_notfull_p7_phi_500
    ),
    title = "500 Hz"
  ),
  list(
    ds = 4,
    pars = list(
      parameters_trials_notfull_p1_phi_250,
      parameters_trials_notfull_p3_phi_250,
      parameters_trials_notfull_p5_phi_250,
      parameters_trials_notfull_p7_phi_250
    ),
    title = "250 Hz"
  )
)

# The test evaluates up to lag 500
lb_lags <- 1:500

ljung_results <- list()

for (trial_n in seq_along(trials)) {
  
  eeg_data <- t(EEGtrial[[trials[trial_n]]])
  cor_mat <- cor(eeg_data)
  
  for (ch in channels_index) {
    
    idx <- order(
      cor_mat[, ch],
      decreasing = TRUE
    )
    
    for (p_index in seq_along(p_values)) {
      
      p <- p_values[p_index]
      
      for (m in models) {

        # downsample function differently applied in function of the dimensions
        if (p == 1) {
          
          y <- eeg_data[, idx[1]]
          y <- downsample(y, m$ds)
          
        } else {
          
          y <- eeg_data[, idx[1:p]]
          
          y <- apply(
            y,
            2,
            downsample,
            m$ds
          )
        }
        
        par <- m$pars[[p_index]][[trial_n]][[ch]]

        # Building function to use in function of the p
        if (p == 1) {
          dlmM1 <- buildSignalnotfullphineq0_p1(par)
        } else if (p == 3) {
          dlmM1 <- buildSignalnotfullphineq0(par)
        } else if (p == 5) {
          dlmM1 <- buildSignalnotfullphineq0(par)
        } else if (p == 7) {
          dlmM1 <- buildSignalnotfullphineq0(par)
        }
        
        eegSmo <- dlmSmooth(y, dlmM1)
        
        signal <- dropFirst(eegSmo$s)
        
        if (p == 1) {
          residuals <- y - signal
        } else {
          residuals <- y[, 1] - signal
        }
        
        # Ljung-Box for lags 1 to 500
        pvals <- sapply(
          lb_lags,
          function(h) {
            Box.test(
              residuals,
              lag = h,
              type = "Ljung-Box"
            )$p.value
          }
        )

        # Store the results
        ljung_results[[length(ljung_results) + 1]] <- data.frame(
          trial = trials[trial_n],
          channel = ordenSensores$electrode[ch],
          p = p,
          fs = m$title,
          lag = lb_lags,
          p_value = pvals,
          significant = pvals < 0.05
        )
      }
    }
  }
}

ljung_df <- bind_rows(ljung_results) # dataframe with 7 columns. The significant column indicates if the test statistic has a p-value lower than 0.05.
