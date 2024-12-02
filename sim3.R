rm(list = ls())
source("share_functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(24)

# Simulation settings
rep_times <- 100
S_sets <- c(2^2, 3^2, 6^2) # Spatial dimension S
TT_sets <- c(20, 50)        # Temporal dimension T
delta_sets <- rbind(c(0,0), c(3,0), c(4,0), c(0,4), 
                    c(0,8), c(3,3), c(4,4)) / 10 # Change size

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

# Define true change point
true_change_point <- T1 # 0.5 * TT

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
  
  # Extract detected change points and p-values
  change_points <- sapply(pettitt_results, function(res) res$change_point)
  change_points <- sapply(change_points, function(x) {
    non_na_values <- x[!is.na(x)]
    if (length(non_na_values) > 0) {
      return(as.numeric(non_na_values[1]))
    } else {
      return(NA)
    }
  })
  
  # Calculate mean, SE, 95% CI, and coverage for detected change points
  mean_cp <- mean(change_points, na.rm = TRUE)
  se_cp <- sd(change_points, na.rm = TRUE) / sqrt(length(change_points))
  ci_lower <- mean_cp - 1.96 * se_cp
  ci_upper <- mean_cp + 1.96 * se_cp
  coverage <- (true_change_point >= ci_lower & true_change_point <= ci_upper)
  
  # Summarize results for this replication
  summary <- data.frame(
    Mean_Change_Point = mean_cp,
    SE = se_cp,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    Coverage = coverage
  )
  
  return(summary)
}

save(final_result, file = "simulation2_results.RData")

# Aggregate results
aggregate_results <- do.call(rbind, final_result)

# Compute overall statistics
overall_coverage <- mean(aggregate_results$Coverage)
overall_mean_change_pt <- mean(aggregate_results$Mean_Change_Point)
overall_se <- mean(aggregate_results$SE)

save(aggregate_results, file = "simulation2_final_results.RData")
