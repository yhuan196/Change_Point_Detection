# Clear workspace and load required libraries
rm(list = ls())
source("share_functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(6)

# Simulation settings
rep_times <- 100
S_sets <- c(2^2, 3^2, 6^2, 8^2) # Spatial dimensions
TT_sets <- c(20, 50, 100, 200)  # Temporal dimensions
delta_sets <- rbind(c(0, 0), c(3, 0), c(4, 0), c(0, 4),
                    c(0, 8), c(3, 3), c(4, 4), c(2, 1)) / 10 # Change sizes

# Function to run simulation for a given setting
run_simulation <- function(S, TT, delta, rep_times) {
  T1 <- T2 <- TT / 2 # Change-point is at 0.5 * TT
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
  results <- foreach(rep_index = 1:rep_times, .packages = c('mvtnorm', 'trend', 'MASS')) %dopar% {
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
    
    # Calculate detection rate
    detection_rate <- mean(!is.na(change_points) & change_points != TT)
    
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
      Coverage = coverage_lambda,
      Detection_Rate = detection_rate
    )
    
    return(summary)
  }
  
  return(results)
}

# Run simulation for both settings
results_list <- list(
  Setting1 = run_simulation(S = S_sets[1], TT = TT_sets[1], delta = delta_sets[8, ], rep_times = rep_times),
  Setting2 = run_simulation(S = S_sets[2], TT = TT_sets[2], delta = delta_sets[8, ], rep_times = rep_times),
  Setting3 = run_simulation(S = S_sets[3], TT = TT_sets[3], delta = delta_sets[8, ], rep_times = rep_times)
)

# Initialize an empty dataframe to store the summary
summary_df <- data.frame(
  T = integer(),
  S = character(),
  Mean_Lambda = numeric(),
  Mean_Coverage = numeric(),
  Mean_CI_Width = numeric(),
  Mean_SE = numeric(),
  Mean_Detection_Rate = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each setting in the results list
for (setting_name in names(results_list)) {
  # Combine the results for the current setting into a single dataframe
  combined_results <- do.call(rbind, results_list[[setting_name]])
  
  # Calculate summary statistics
  mean_lambda <- mean(combined_results$Estimated_Lambda, na.rm = TRUE)
  mean_coverage <- mean(combined_results$Coverage, na.rm = TRUE)*100
  mean_ci_width <- mean(combined_results$CI_Upper - combined_results$CI_Lower, na.rm = TRUE)
  mean_ci_lower <- mean(combined_results$CI_Lower, na.rm = TRUE)
  mean_ci_upper <- mean(combined_results$CI_Upper, na.rm = TRUE)
  mean_se <- mean(combined_results$SE_Lambda, na.rm = TRUE)
  mean_detection_rate <- mean(combined_results$Detection_Rate, na.rm = TRUE)*100
  
  # Add the results to the summary dataframe with math symbols in S
  summary_df <- rbind(summary_df, data.frame(
    T = if (setting_name == "Setting1") {
      20
    } else if (setting_name == "Setting2") {
      50
    } else { # Setting3
      100
    },
    S = if (setting_name == "Setting1") {
      "2^2"
    } else if (setting_name == "Setting2") {
      "3^2"
    } else { # Setting3
      "2^3"
    }, # Math symbols as strings
    Mean_Lambda = mean_lambda,
    Mean_Detect_Perc = mean_detection_rate,
    Mean_Coverage_Perc = mean_coverage,
    Mean_SE = mean_se,
    Mean_CI_Lower = mean_ci_lower,
    Mean_CI_Upper = mean_ci_upper,
    Mean_CI_Width = mean_ci_width
    
  ))
  
}

# Print the summary dataframe
print(summary_df)
save(summary_df, file = "pettitt_table2.RData")


