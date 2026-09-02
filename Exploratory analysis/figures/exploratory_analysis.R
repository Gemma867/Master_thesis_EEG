
# This .R file contains the code to plot the figures of the exploratory analysis.

library(dplyr)
library(ggplot2)
library(tidyr)
library(zoo)
library(readxl)

# sensor names and order corresponding to the order of EEG data list
ordenSensores <- read_excel("ordenSensores.xls", col_names = FALSE)
names(ordenSensores)<-c("electrode")
# EEG data list
EEGtrial <- readRDS("matrices.rds")

# 2d coordinates of the electrodes locations (extracted from the electrode layout scheme)
chan_coords <- data.frame(
  electrode = ordenSensores$electrode,
  x = c(0.42667,0.67556,0.65778,0.46222,0.19556,0.16889,0.56889,0.86222,0.99556,0.92444,
        0.70222,0.46222,0.19556,-0.07111,-0.15111,-0.05333,0.24889,0.58667,0.96,1.24444,
        1.34222,1.12,0.80889,0.45333,0.07111,-0.24,-0.5,-0.43556,-0.15111,0.17778,
        0.41778,1.02222,1.42222,1.64444,1.51111,1.06667,0.65778,0.22222,-0.17778,-0.64889,
        -0.80889,-0.57778,-0.20444,1.81333,2,1.82222,1.39556,0.93333,0.43556,-0.00889,
        -0.5,-0.95111,-1.18222,-0.96,2.24889,2.00889,1.42222,0.8,0.17778,-0.51556,-1.15556,-1.5,-0.39111,1.44),
  y = c(0.48,0.32889,0.00889,-0.12444,-0.00889,0.28444,0.78222,0.58667,0.19556,-0.09778,
        -0.31111,-0.4,-0.32,-0.09778,0.21333,0.56889,0.77333,1.07556,0.92444,0.54222,
        0.05333,-0.39111,-0.60444,-0.68444,-0.63111,-0.39111,0,0.50667,0.92444,1.07556,
        1.39556,1.25333,0.93333,0.31111,-0.41778,-0.77333,-0.96889,-1.00444,-0.82667,-0.41778,
        0.29333,0.94222,1.23556,1.05778,0.5,-0.62222,-0.96,-1.18222,-1.25333,-1.22667,
        -1,-0.67556,0.47111,0.93333,0.68444,-0.88,-1.31556,-1.55556,-1.57333,-1.36,-0.84444,0.5,1.84889,1.71556)
)

#### Boxplot of mean general statistics per channel (1st quartile, 3rd quartile, minimum, maximum, median)

def_total_real_list <- readRDS("def_total_real_list.rds") # List of 64 elements (channels) with dataframes of 3000 rows (samples) and 761 columns (time + 760 trials). Obtained from EEGtrial.
general_statistics <-data.frame(channel=ordenSensores$electrode) # list to store the mean IQR

for (i in 1:64){
  q25<-mean(apply(def_total_real_list[[i]][, -1], 2, function(x) quantile(x, 0.25)))
  q75<-mean(apply(def_total_real_list[[i]][, -1], 2, function(x) quantile(x, 0.75)))
  min<-mean(apply(def_total_real_list[[i]][, -1], 2, min))
  max<-mean(apply(def_total_real_list[[i]][, -1], 2, max))
  median<-mean(apply(def_total_real_list[[i]][, -1], 2, median))
  IQR<-q75-q25
  
  general_statistics$IQR[i]<-IQR
  general_statistics$q25[i] <- q25
  general_statistics$q75[i] <- q75
  general_statistics$min[i] <- min
  general_statistics$max[i] <- max
  general_statistics$median[i] <- median
}

ggplot(general_statistics, aes(x = factor(channel, levels=channel))) +
  geom_boxplot(
    aes(
      ymin = min,
      lower = q25,
      middle = median,
      upper = q75,
      ymax = max
    ),
    stat = "identity",
    fill = "mediumturquoise",
    color = "darkblue"
  ) +
  xlab("Channel") +
  ylab("Amplitude") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(
    floor(min(general_statistics$min, na.rm = TRUE)),
    ceiling(max(general_statistics$max, na.rm = TRUE)),
    by = 5
  ))

#### Boxplot of standard deviations per channel

# Standard deviations computation for every trial and channel
var_df <- lapply(EEGtrial, function(mat) apply(mat, 1, sd))     # compute row means for each matrix
var_df <- as.data.frame(var_df)

rownames(var_df)<-ordenSensores$electrode

df_t <- as.data.frame(t(var_df))
df_t$original_col <- rownames(df_t)

# Convert dataframe to plot
df_long_var <- pivot_longer(
  df_t,
  cols = -original_col,
  names_to = "channel",
  values_to = "value"
)

# Compute medians
medians <- df_long_var %>%
  group_by(channel) %>%
  summarise(med = median(value), .groups = "drop")
# Join medians
df_long_var <- df_long_var %>%
  left_join(medians, by = "channel")

# Reorder channel factor by median
df_long_var$channel <- reorder(df_long_var$channel, df_long_var$med)

ggplot(df_long_var, aes(x = channel, y = value, fill = med)) +
  geom_boxplot(width = 0.75) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.75,
    color = "red",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  labs(
    x = "Electrode",     
    y = "Standard deviation",       
    fill = "Median std"   
  ) +
  theme_bw()+
  scale_fill_gradient2(low = "blue",mid="lightgreen",high = "red",midpoint=4.5,    limits = c(2, 8.5), oob = scales::squish  ) +
  coord_cartesian(ylim = c(1, 10)) +
  scale_y_continuous(breaks = seq(1, 10, by = 1)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size=11),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size = 15),  # bigger x label
    axis.title.y = element_text(size = 15)   # bigger y label
  )

#### Electrode layout median variance plot

library(ggrepel)

chan_coords$med_var <- medians$med[match(chan_coords$electrode, medians$channel)]

ggplot(chan_coords, aes(x = x, y = y, color = med_var)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = electrode), size = 3) + 
  scale_color_gradient2(
    low = "blue",
    mid = "lightgreen",
    high = "red",
    midpoint = 4.5,
    limits = c(2, 8.5),
    oob = scales::squish
  ) +
  coord_cartesian(ylim = c(1, 10)) +
  labs(color = "Median std") +
  coord_fixed() +
  theme_void()

#### Some EEG signals examples

channel_indices <- c(14,1,44,34,54,42,61,63,64)
trial<-39
df <- EEGtrial[[trial]][channel_indices, ]
df$electrode <- factor(
  ordenSensores$electrode[channel_indices],
  levels = ordenSensores$electrode[channel_indices]
)

df$channel <- 1:nrow(df)

# Convert to long format
df_long <- df %>%
  pivot_longer(-c(channel, electrode), names_to = "time", values_to = "amplitude")

df_long$time <- as.numeric(gsub("V", "", df_long$time))

# Add vertical offset 
df_long <- df_long %>%
  group_by(channel) %>%
  mutate(amplitude_offset = amplitude)

# Set y-axis limits
my_ymin <- -20
my_ymax <- 20

ggplot(df_long, aes(x = time, y = amplitude_offset)) +
  geom_line(color = "blue") +
  facet_wrap(~electrode, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_continuous(limits = c(0, 3000)) +
  #scale_y_continuous(limits = c(my_ymin, my_ymax)) +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "green", linewidth = 0.8) +
  theme_minimal() +
  labs(
    title = paste0("EEG Channels trial ", trial),
    x = "Time (ms)",
    y = "Amplitude"
  ) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  )

#### Example of electrode layout correlation plot

cor_df <- lapply(matrices, function(mat) cor(t(mat)))
# Fisher Z transformation
z_df <- lapply(cor_df, atanh)
# Mean Fisher Z
mean_z <- Reduce("+", z_df) / length(z_df)
# Transform back to correlations
mean_cor <- tanh(mean_z)

chan_coords$corFP1 <- mean_cor[43,]

# Fp1 mean correlation values with other electrodes plot
ggplot(chan_coords, aes(x = x, y = y, color = corFP1)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = electrode), size = 3) + 
  scale_color_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint=0,
    oob = scales::squish
  )+
  labs(color = "Average correlation",
       title = "Fp1")+
  coord_fixed() +
  theme_void()

