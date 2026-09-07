
# First model:

buildSignal <- function(para) {
  
  # Construct the lower-triangular matrix L from the first
  # (p^2 + p)/2 elements of the parameter vector.
  # The triangular structure ensures that V = L'L is positive semidefinite.
  L <- matrix(0, p, p)
  L[upper.tri(L, TRUE)] <- para[1 : ((p*p+p)/2)]
  
  # Exponentiate the diagonal elements of L to ensure that
  # the resulting covariance matrix is positive definite.
  diag(L) <- exp(diag(L))
  
  modSignal <- dlm(
    FF = matrix(1, p, 1), # Observation matrix
    V = crossprod(L), # Observation noise covariance matrix
    GG = 1, # State transition matrix (local level)
    W = exp(para[(p*p+p)/2+1]), # State noise variance
    m0 = 0, # Initial state mean
    C0 = 1e7 # Initial state variance set to a very large value
  )
  
  return(modSignal)
}

# Second model:

# The definition of the model follows the same order.
  
# p > 1:

buildSignalnotfullphineq0 <- function(para) {
  
  # Initialize the observation noise covariance matrix V.
  # The first channel is treated separately, with its variance
  # parameterized independently.
  VV <- matrix(0, p, p)
  VV[1, 1] <- exp(para[1])

  # Construct the covariance matrix for the remaining p-1 channels.
  # The upper-triangular elements of L are estimated from the
  # parameter vector. The diagonal elements are exponentiated
  # to ensure positive diagonal entries.
  L <- matrix(0, p-1, p-1)
  L[upper.tri(L, TRUE)] <- para[2 : (1+((p-1)*(p-1)+(p-1))/2)]
  diag(L) <- exp(diag(L))
  LL <- crossprod(L) 
  
  # Obtain the covariance matrix.
  VV[2:p, 2:p] <- LL
  
  ind <- (1+((p-1)*(p-1)+(p-1))/2)
  
  # Construct the observation matrix F.
  # The first channel has a fixed loading of 1, while the
  # remaining p-1 loadings are estimated.
  FF <- matrix(1, p, 1)
  FF[2:p,1] <- para[(ind+1):(ind+p-1)]
  
  ind <- (ind+p-1)
  
  modSignal <- dlm(
    FF = FF, 
    V = VV+matrix(1,p,p)*1e-6,
    GG = exp(para[ind+2]), 
    W = exp(para[ind+1]), 
    m0 = 0, 
    C0 = 1e4
  )
  
  return(modSignal)
}

# p = 1:

buildSignalnotfullphineq0_p1 <- function(para) {
  
  VV <- matrix(0, p, p)
  VV[1, 1] <- exp(para[1])
  
  FF <- 1
  
  modSignal <- dlm(
    FF = FF, 
    V = VV+1e-6,
    GG = para[2], 
    W = exp(para[3]), 
    m0 = 0, 
    C0 = 1e4
  )
  
  return(modSignal)
}