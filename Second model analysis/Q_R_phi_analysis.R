
# This .R file contains the code to obtain the boxplots of phi, Q / R, Var(x) / R, and the scatter plot of R_250/R_500 against Q_250/Q_500.

library(ggplot2)
library(dplyr)
library(gsignal)
library(dlm)
library(tidyr)

parameters_trials_notfull_p1_phi<-readRDS("parameters_trials_notfull_p1_phi.rds")
parameters_trials_notfull_p1_phi_250<-readRDS("parameters_trials_notfull_p1_phi_250.rds")
parameters_trials_notfull_p1_phi_500<-readRDS("parameters_trials_notfull_p1_phi_500.rds")
source("Functions.R")
EEGtrial <- readRDS("matrices.rds")

p<-1

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

ratio_df <- data.frame()

for (trial_n in 1:6){
  
  eeg_data <- t(EEGtrial[[trials[trial_n]]])
  
  for (ch in channels_index) {
    
    idx <- order(cor(eeg_data)[, ch], decreasing = TRUE)
    
    for (m in models) {
      
      para<-m$pars[[trial_n]][[ch]]
      
      Q<-exp(para[3])
      Q_r <- exp(para[3])/(1-para[2]^2)
      R<-exp(para[1])
      phi<-para[2]
      ratio<-Q/R
      ratio_r <- Q_r/R
      ratio_df <- rbind(
        ratio_df,
        data.frame(
          channel = ch,
          ratio = ratio,
          ratio_r=ratio_r,
          phi = phi,
          Q=Q,
          R=R,
          model = m$title,
          trial = trials[trial_n]
        ))
      
      
    }
  }
}

ratio_df$model <- factor(
  ratio_df$model,
  levels = c("1000 Hz", "500 Hz", "250 Hz")
)

#### Scatter plot R_250/R_500 against Q_250/Q_500

ratio_plot <- ratio_df %>%
  dplyr::filter(model %in% c("500 Hz", "250 Hz")) %>%
  dplyr::filter(ratio < 30) %>%
  dplyr::select(trial, Q, R, model, channel) %>%
  tidyr::pivot_wider(
    names_from = model,
    values_from = c(Q, R)
  ) %>%
  dplyr::mutate(
    Q_ratio = `Q_250 Hz` / `Q_500 Hz`,
    R_ratio = `R_250 Hz` / `R_500 Hz`
  )

ggplot(ratio_plot, aes(x = Q_ratio, y = R_ratio)) +
  geom_point() +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  labs(
    x = expression(Q[250~Hz] / Q[500~Hz]),
    y = expression(R[250~Hz] / R[500~Hz])
  ) +
  theme_minimal()

### Boxplots

# Filter out two outliers (ratio = 32 and ratio = 299)
ratio_plot_df <- ratio_df %>%
  dplyr::filter(ratio < 30) %>%
  dplyr::group_by(model) %>%
  dplyr::mutate(
    min_ratio = min(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE)
  )

# Boxplot phi

ggplot(ratio_plot_df, aes(x = model, y = phi, fill = model)) +
  geom_boxplot() +
  labs(
    x = "Sampling frequency",
    y = "Q / R",
    fill = "fs"
  ) +
  theme_minimal()+
  theme(
    legend.position = "none"
  )

# Boxplot process variance / R

ggplot(ratio_plot_df, aes(x = model, y = ratio_r, fill = model)) +
  geom_boxplot() +
  labs(
    x = "Sampling frequency",
    y = "Variance of process / R",
    fill = "fs"
  ) +
  theme_minimal()+
  theme(
    legend.position = "none"
  )

# Boxplot Q / R

ggplot(ratio_plot_df, aes(x = model, y = ratio, fill = model)) +
  geom_boxplot() +
  labs(
    x = "Sampling frequency",
    y = "Q / R",
    fill = "fs"
  ) +
  theme_minimal()+
  theme(
    legend.position = "none"
  )

