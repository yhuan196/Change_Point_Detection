rm(list = ls())
source("share_functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(24)

# Simulation settings
rep_times <- 100
S_sets <- c(2^2, 3^2, 6^2) # Spatial dimensions
TT_sets <- c(20, 50)        # Temporal dimensions
delta_sets <- rbind(c(0,0), c(3,0), c(4,0), c(0,4), 
                    c(0,8), c(3,3), c(4,4)) / 10 # Change sizes

# Select a specific setting
S <- S_sets[1]
TT <- TT_sets[2]
T1 <- T2 <- TT / 2 # Change-point is at 0.5 * TT
delta <- delta_sets[3, ]

# True model parameters
theta1 <- c(-0.4, 0.5, 1)
theta2 <- c(-0.2 + delta[1], 0.5 + delta[2], 1)

# Generate spatial grid and distance matrix
min.dist <- 1
max.lati <- max.long <- sqrt(S)
mat.lattice <- cbind(rep(seq(1, max.lati, by = min.dist), each = sqrt(S)),
                     rep(seq(1, max.long, by = min.dist), sqrt(S)))
S.dist <- round(as.matrix(dist(mat.lattice, upper = TRUE, diag = TRUE)), 10)

# Define true lambda
true_lambda <- T1 / TT

# Main simulation loop
final_result <- foreach(rep_index = 1:rep_times, .packages = c('mvtnorm', 'trend', 'MASS')) %dopar% {
  set.seed(1128 + rep_index)
  # Simulate dataset with one change-point
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
  
  # Extract detected change points
  change_points <- sapply(pettitt_results, function(res) res$change_point)
  change_points <- sapply(change_points, function(x) {
    non_na_values <- x[!is.na(x)]
    if (length(non_na_values) > 0) {
      return(as.numeric(non_na_values[1]))
    } else {
      return(NA)
    }
  })
  
  # Calculate lambda as the average detected change point divided by TT
  estimated_change_point <- mean(change_points, na.rm = TRUE)
  estimated_lambda <- estimated_change_point / TT
  
  se_lambda <- sd(change_points, na.rm = TRUE) / (sqrt(length(change_points)) * TT)
  ci_lower_lambda <- estimated_lambda - 1.96 * se_lambda
  ci_upper_lambda <- estimated_lambda + 1.96 * se_lambda
  coverage_lambda <- (true_lambda >= ci_lower_lambda & true_lambda <= ci_upper_lambda)
  
  summary <- data.frame(
    Estimated_Lambda = estimated_lambda,
    SE_Lambda = se_lambda,
    CI_Lower = ci_lower_lambda,
    CI_Upper = ci_upper_lambda,
    Coverage = coverage_lambda
  )
  
  return(summary)
}

# Save results
save(final_result, file = "simulation_lambda_results.RData")

aggregate_results <- do.call(rbind, final_result)

mean(aggregate_results$Estimated_Lambda)
mean(aggregate_results$SE_Lambda)
mean(aggregate_results$Coverage)
