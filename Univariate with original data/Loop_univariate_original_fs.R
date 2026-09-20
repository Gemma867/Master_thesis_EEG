
# This .R file contains the code to obtain the parameter estimation of the univariate models with sampling frequency of 1000 Hz.

library(GGally)
library(ggplot2)
library(tidyr)
library(zoo)
library(readxl)
library(dplyr)
library(gsignal)
library(dlm)

# Local level:

source("Functions_univariate_original.R")
EEGtrial<-readRDS("matrices.rds")

# Estimation of parameters for 760 trials and 64 channels

params_list_def <- vector("list", 64) # list with parameters. Each element corresponds to the estimation of the 760 trials for each channel.

# Loop for the 64 channels
for (j in 1:64) {
  params_j <- matrix(NA, nrow = 760, ncol = 2)
  
  # Loop for the 760 trials
  for (i in 1:760) {
    matrix <- EEGtrial[[i]]
    EEG <- matrix[j, ]
    
    # Maximum likelihood estimation of the parameters
    fit <- dlmMLE(as.numeric(EEG), parm = c(0,0), build = buildFun1)
    
    params_j[i, ] <- fit$par
    
  }
  
  # Store the parameters
  params_list_def[[j]] <- params_j
}

# Univariate AR(1) process:

params_univariate_process_matrix <- vector("list", 64)

for (j in 1:64) {
  params_j <- matrix(NA, nrow = 760, ncol = 2)
  
  # Loop for the 760 trials
  for (i in 1:760) {
    
    matrix <- EEGtrial[[i]]
    EEG <- matrix[j, ]
    
    fit <- dlmMLE(as.numeric(EEG), parm = c(-10,-10,1), build = buildFun2)
    
    params_j[i, ] <- fit$par
    
  }
  
  # Store the parameters
  params_univariate_process_matrix[[j]] <- params_j
}

