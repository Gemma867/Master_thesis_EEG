
# This .R file contains the code to the functions to build the univariate filters.

library(dlm)

# Local level

buildFun1 <- function(x) {
  dlmModPoly(1, dV = exp(x[1]), 
             dW = exp(x[2]))} # only measurement and process variances estimated


# Univariate AR(1) process

buildFun2 <- function(x) {
  dlm(
    FF = matrix(1),          # Observation matrix (F)
    V  = exp(x[1]),          # Observation variance
    GG = matrix(x[3]),          # Process matrix (G)
    W  = exp(x[2]),          # Process variance
    m0 = 0,
    C0 = 1e7
  )
}