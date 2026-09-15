
# This .R file contains the code to obtain the plots from the unvivariate results analysis with the original frequency and over all the trials (760) and channels (64) available.

library(ggplot2)
library(dplyr)

# Boxplot of measurement or process noise covariance per channel

# Create long-format data frame from the estimated parameters
param_list_long <- lapply(seq_along(params_univariate_process_matrix_def), function(i) {
  
  mat <- params_univariate_process_matrix_def[[i]]
  
  data.frame(
    value = c(
      exp(mat[, 1]),  # Measurement noise covariance R
      exp(mat[, 2]),  # Process noise covariance Q
      mat[, 3]        # AR(1) parameter phi
    ),
    
    column = factor(
      rep(
        c("Measurement", "Process", "Process matrix"),
        each = nrow(mat)
      ),
      levels = c("Measurement", "Process", "Process matrix")
    ),
    
    electrode = factor(i),
    trial = rep(1:nrow(mat), times = 3)
  )
})

# Combine all electrodes. This dataframe will be used for the following plots also.
df_param_list_long <- do.call(rbind, param_list_long)

# Divide electrodes into two panels
df_param_list_long$panel <- ifelse(
  as.numeric(df_param_list_long$electrode) <= 32,
  "1-32",
  "33-64"
)

# Add electrode names
df_param_list_long$electrode_name <- ordenSensores$electrode[
  as.numeric(df_param_list_long$electrode)
]

# Preserve the electrode order
df_param_list_long$electrode_name <- factor(
  df_param_list_long$electrode_name,
  levels = ordenSensores$electrode
)

p1 <- ggplot(
  df_param_list_long[
    #df_param_list_long$value <= 100 & # Filter the dataframe if you want to take out outliers
      df_param_list_long$column == "Measurement", # Choose Measurement or Process in function of the plot wanted
  ],
  aes(x = electrode_name, y = value)
) +
  geom_boxplot(
    outlier.colour = "red2",
    outlier.shape = 16,
    outlier.size = 1.5
  ) +
  facet_wrap(~ panel, scales = "free_x", ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(
      size = 14,
      face = "bold",
      margin = margin(r = 10)
    ),
    strip.text = element_blank()
  ) +
  labs(
    x = "Electrodes",
    y = "Measurement noise covariance estimate"
  )

p1

# Barplot of medians per channel

library(dplyr)

# Medians of R stored in a datframe
df_median_proc <- df_param_list_long %>%
  dplyr::filter(column == "Process") %>%
  dplyr::group_by(panel, electrode_name) %>%
  dplyr::summarise(median_value = median(value), .groups = "drop")

# Barplot of the median values
p3 <- ggplot(df_median_proc, aes(x = electrode_name, y = median_value)) +
  geom_col(fill = "violetred") +   
  facet_wrap(~ panel, scales = "free_x", ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    strip.text = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  #ylim(0,10) +
  labs(
    x = "Electrodes",
    y = "Median Process noise covariance estimate"
  )

p3

# Boxplot of process matrix estimate per channel

p4<-ggplot(
  df_param_list_long[(df_param_list_long$column == "Process matrix"), ],
  aes(x = electrode_name, y = value)
) +
  geom_boxplot(
    fill = "skyblue",
    colour = "darkblue",
    outlier.colour = "darkblue",      # lilac/purple color for outliers
    outlier.shape = 16,           # solid circle
    outlier.size = 1.5  
  ) +
  facet_wrap(~ panel, scales = "free_x", ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    strip.text = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Electrodes", y = "Process matrix estimate", fill = NULL)

p4


# Variance reduction boxplot

library(dlm)

# Store variance reduction for each electrode
var_reduction_uni <- vector("list", 64)

for (j in 1:64) {
  
  # Matrix to store VRR for the 760 trials
  var_reduction_j <- numeric(760)
  
  for (i in 1:760) {
    
    # Original EEG signal
    EEG <- as.numeric(EEGtrial[[i]][j, ])
    
    # Estimated parameters for this electrode and trial. Change the list if you want the local level model.
    para <- params_univariate_process_matrix_def[[j]][i, ] 
    
    # Build the state-space model using the already estimated parameters
    model <- buildFun2(para)
    
    # Kalman smoother
    smoothed <- dlmSmooth(EEG, model)
    
    # Extract the smoothed latent signal
    denoised <- dropFirst(smoothed$s)
    
    # Variance reduction
    var_reduction_j[i] <- (
      var(EEG) - var(denoised)
    ) / var(EEG)
  }
  
  # Store results for this electrode
  var_reduction_uni[[j]] <- var_reduction_j
}



var_reduction_long <- do.call(rbind, lapply(seq_along(var_reduction_uni), function(i) {
  data.frame(
    value = var_reduction[[i]],       
    group = ordenSensores$electrode[i]
  )
}))

first_half <- ordenSensores$electrode[1:32] # dividing for plot clarity
second_half <- ordenSensores$electrode[33:64]

var_reduction_long$panel <- ifelse(
  var_reduction_long$group %in% first_half,
  "First 32",
  "Last 32"
)

var_reduction_long$group <- factor(
  var_reduction_long$group,
  levels = ordenSensores$electrode
)

# Plot with y-axis 0 to 100%
p_var_1 <- ggplot(var_reduction_long, aes(x = group, y = value)) +
  geom_boxplot(outlier.colour = "darkorange1") +
  facet_wrap(~panel, scales = "free_x", ncol = 1) +
  theme_minimal()+
  labs(
    x = "Electrode",
    y = "Variance reduction"
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)))

# Plot with zoom in y-axis
p_var_2 <- ggplot(var_reduction_long[var_reduction_long$value < 0.5,], aes(x = group, y = value)) +
  geom_boxplot(outlier.colour = "darkorange1") +
  facet_wrap(~panel, scales = "free_x", ncol = 1) +
  theme_minimal()+
  labs(
    x = "Electrode",
    y = "Variance reduction"
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)))

p_var_1

p_var_2
# 