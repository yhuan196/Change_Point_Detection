rm(list=ls())
source("share_functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(24)
############# Code for simulation 1


### All simulation settings considered in Simulation 1
rep_times <- 100
S_sets <- c(2^2, 3^2, 6^2) # Spatial dimension S
TT_sets <- c(20, 50) # Temporal dimension T
delta_sets <- rbind(c(0,0), c(3,0), c(4,0), c(0,4), 
                    c(0,8), c(3,3), c(4,4))/10 # Change size

### Select a specific setting in Simulation 1 to run
# i.e. combination of (S, T, delta)
S <- S_sets[1]
TT <- TT_sets[1]
T1 <- T2 <- TT/2 # segment length (change-point is at 0.5*TT)
delta <- delta_sets[3,]

# True model parameters in the first and second segment c(phi, rho, sigma2)
theta1 <- c(-0.4, 0.5, 1)
theta2 <- c(-0.2 + delta[1], 0.5 + delta[2], 1)

### Generate matrix of 2D regular grid and distance matrix
min.dist <- 1
max.lati <- max.long <- sqrt(S)
n.lati <- max.lati/min.dist
n.long <- max.long/min.dist
mat.lattice <- cbind(rep(seq(1, max.lati, by=min.dist), n.long), 
                     rep(seq(1, max.long, by=min.dist), each=n.lati))
# round to 10 decimal place to remove extra distance combinations
S.dist <- round(as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE)), 10) 

### Parameters for implementing CLDML
# Set spatial and temporal lag in computing CL
s.lag <- t.lag <- 1

######### Main function ##########
# Main simulation loop
final_result <- foreach(rep_index = 1:rep_times, .packages = c('mvtnorm','MASS', 'trend')) %dopar% {
  # Simulate a dataset with one single change-point at 0.5 * TT
  y <- rbind(sim.y(theta = theta1, S.dist = S.dist, TT = T1, T.burn = 100),
             sim.y(theta = theta2, S.dist = S.dist, TT = T2, T.burn = 100))
  
  # Apply Pettitt's test to each spatial location
  pettitt_results <- apply(y, 2, function(ts) {
    result <- pettitt.test(ts)
    list(
      change_point = result$estimate,
      p_value = result$p.value
    )
  })
  
  # Extract change points and p-values
  change_points <- sapply(pettitt_results, function(res) res$change_point)
  change_points <- sapply(change_points, function(x) {
    non_na_values <- x[!is.na(x)]
    if (length(non_na_values) > 0) {
      return(as.numeric(non_na_values[1]))
    } else {
      return(NA)
    }
  })
  p_values <- sapply(pettitt_results, function(res) res$p_value)
  
  # Summarize results for this replication
  summary <- data.frame(
    Location = 1:ncol(y),
    Change_Point = change_points,
    P_Value = p_values
  )
  
  return(summary)
}


save.image("Simulation1_.RData")
final_result
save(final_result, file = "Simulation1_Pettitt1.RData")


