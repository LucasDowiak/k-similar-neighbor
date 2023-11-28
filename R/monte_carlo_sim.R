setwd("~/Git/dtw-in-finance")
source("R/similarity_functions.R")
library(data.table)
library(mvtnorm)

B <- 500
N <- 231
var1 <- 0.0009252
var2 <- 0.0009252
rho <- 0.25
nu <- 10


# For given parameter values simulate multivariate
simulate_cor_vs_dtw <- function(B, N, var1, var2, rho, nu)
{
  mu <- c(0, 0)
  covariance <- sqrt(var1) * sqrt(var2) * rho
  sigma <- matrix(c(var1, rep(covariance, 2), var2), nrow=2)
  pairs <- list(list("V1", "V2")); names(pairs) <- "V1"
  
  # Sample returns from the multivariate t-distribution
  samples <- replicate(B, rmvt(n=N - 1L, sigma=sigma, df=nu, delta=mu), simplify=FALSE)
  
  # Function to transform returns into their standardized levels: y_t = cumprod(exp(r_t))
  prep_data <- function(M)
  {
    out <- apply(M, MARGIN=2, FUN=function(x) c(1, cumprod(exp(x))))
    return(as.data.table(out))
  }
  samples <- lapply(samples, prep_data)
  
  # Run dtw 
  results <- lapply(samples, calculate_similarity, metric=dynamic_time_warp, pairs=pairs, ncores=1L)
  
  # Return data.table of results
  dtws <- sapply(results, function(x) as.numeric(x[[1]]))
  return(data.table(B=B, N=N, var1=var1, var2=var2, rho=rho, nu=nu, dtw=dtws))
}

tst <- simulate_cor_vs_dtw(B, N, var1, var2, rho, nu)
system.time(simulate_cor_vs_dtw(B, N, var1, var2, rho, nu))

RHO <- seq(-0.4, 1, 0.01)
VAR <- seq(0.0001, 0.012, 0.0001)
NU <- seq(2, 15, 1)
length(VAR)**2 * length(RHO) * length(NU)
