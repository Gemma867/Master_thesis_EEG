
# This .R file contains the code to estimate the parameters of the Kalman filter with function buildSignal() for the 6 trials selected and the 64 channels.

library(dlm)
source("Functions.R")

trials <- c(1,2,4,5,21,223)
p <- 7
parameters_6trials_junt <- list()
EEGtrial <- readRDS("matrices.rds")

npar <-(p*p + p)/2 + 1

for (trial in trials){
  
  parameters_6trials_junt[[as.character(trial)]] <- list()
  
  eeg_data <- t(EEGtrial[[trial]]) # select the EEG data from trial
  
  for (ch in 1:64){
    
    cat("\nTrial:", trial, "Channel:", ch, "-", chan_names[ch])
    
    idx <- order(cor(eeg_data)[,ch], decreasing = TRUE) # extract index of the channels ordered in function of correlation with reference channel
    
    for (i in 1:p){
      cat("\n    Channel:", idx[i])
    }
    
    y <- eeg_data[, idx[1:p]] # Select p signals
    colnames(y) <- chan_names[idx[1:p]]
    
    # Maximum Likelihood estimation
    fitSignal <- dlmMLE(
      y,
      parm = rep(0, npar),
      build = buildSignal,
      hessian = TRUE,
      control = list(maxit = 500)
    )
    
    parameters_6trials_junt[[as.character(trial)]][[ch]] <- fitSignal$par
  }
}