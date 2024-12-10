#'---------------------------------------------------------------------------------------------------------
#' Project Name: utilities_main.R
#' 
#' Short Description: Run simution for petitt's test and CLMDL
#'
#' Author: Yi Huang, James Robertson 
#'
#' Data Authored: 11/28/2024
#'---------------------------------------------------------------------------------------------------------
#'

# Clear workspace and load required libraries
rm(list = ls())
source("utilities_generate.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(6)

#'---------------------------------------------------------------------------------------------------------
#' Section: Simulation for Petitt Test
#' PURPOSE: 
#' @param rep_times
#' @param S_sets
#' @param TT_sets
#' @param delta_set
#'---------------------------------------------------------------------------------------------------------
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





#'---------------------------------------------------------------------------------------------------------
#' Section: Simulation for CLMDL
#' PURPOSE: 
#' @param Q
#' @param y
#' @param S.dist
#' @param S
#'---------------------------------------------------------------------------------------------------------
# ### Read-in the simulated spatio-temporal process
# # data_name_set <- c('y0','y1','y2')
# # data_name <- 'y0' # no change-point
# data_name <- 'y1' # one change-point
# # data_name <- 'y2' # two change-point
# 
# 
# data <- readRDS(paste0(data_name, '.RDS'))
# y <- data$y # T times S spatio-temporal data
# S.dist <- data$S.dist # distance matrix of the S spatial locations
source("utilities_estimate.R")
Q<-500
for (siz in 'large') {#c('small','med','large')) {
  timStart<-proc.time()
  savResults<-foreach(q=1:Q,
                    .inorder=F,
                    .packages = c('mvtnorm')) %dopar%  {
                      print(siz)
print(paste('iter',q))
set.seed(1128 + q)

# set.seed(500)
smalSamp<-T
source('spatial_ar1_data.R',local=T)
y<-result$data
S.dist<-result$distance_matrix
S <- ncol(y) # spatial dimension
TT <- nrow(y) # time dimension

### Parameters for implementing CLDML
# Set spatial and temporal lag in computing CL
s.lag <- t.lag <- 1

# Collect all pairs of observations within s.lag distance (will be used in PL)
s.comb <- s.dist.pair(S.dist, s.lag) # time.lag > 0
s.comb.0 <- s.dist.pair.0(S.dist, s.lag) # time.lag = 0
s.len <- unique(s.comb[,3])
s.len <- s.len[s.len>0] # all unique spatial distance

# Calculate average number of times an obs used in CL for Ck
pair_stat <- D.cal(y, S.dist, s.lag, t.lag)
Ck <- mean(pair_stat$use.obs)
remedy <- pair_stat$remedy # number of marginal likelihood needed for correcting edge effect

# initial parameter estimation based on entire data using the 4-parameter spatial auto-regressive model
ini <- optim(c(0,1,1,0), pl, y=y, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
             method="L-BFGS-B", lower=c(-0.7,0.1,1e-3,-1), upper=c(0.7,3,5,1), ts=F)$par

# Calculate K for pruning step
p.min <- p.max <- length(ini)
K <- floor(Ck)*(floor(log(S*TT))*(p.min/2 - p.max) + (2 + p.max)*log(2) - floor(log(TT)))

######### Main function ##########
### Run PELT for minimizing CLMDL (this is the default algorithm in the paper.)
t0 <- proc.time()
pelt_result <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K,
                    t.lag=t.lag, min.length=0.1*TT, ini=ini, res=1)
t_compute <- proc.time() - t0

### One can run PELT with res=2 to further speed up computation
# (res=2 means we assume that cp can only happen on t such that t%%2==0, i.e. t is even.)
# Therefore the resolution of this algorithm is 2.
# t0 <- proc.time()
# pelt_result2 <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K,
#                      t.lag=t.lag, min.length=0.1*TT, ini=ini, res=2)
# t_compute2 <- proc.time() - t0

### Subsequent analysis after CP detection
num_cp <- length(pelt_result)-2
ci_result <- c()
if(num_cp>0){
  for(cp_index in 1:num_cp){
    # Parameter estimation on the pre-change segment
    y.seg1 <- y[(pelt_result[cp_index]+1):pelt_result[cp_index+1],]
    pmle1 <- optim(ini, pl, y=y.seg1, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
                   method="L-BFGS-B", lower=c(-0.7,0.1,1e-3,-1), upper=c(0.7,3,5,1), ts=F)
    # Parameter estimation on the post-change segment
    y.seg2 <- y[(pelt_result[cp_index+1]+1):pelt_result[cp_index+2],]
    pmle2 <- optim(ini, pl, y=y.seg2, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
                   method="L-BFGS-B", lower=c(-0.7,0.1,1e-3,-1), upper=c(0.7,3,5,1), ts=F)
    
    seg1_tmp <- nrow(y.seg1)
    seg2_tmp <- nrow(y.seg2)
    TT_tmp <- seg1_tmp + seg2_tmp
    
    # CI construction
    B <- 100
    theta1_est <- pmle1$par
    theta2_est <- pmle2$par
    q_emp <- c()
    t0 <- proc.time()
    for(b in 1:B){
      set.seed(b)
      print(b)
      yB <- rbind(sim.y(theta=theta1_est, S.dist=S.dist, TT=seg1_tmp, T.burn=100),
                  sim.y(theta=theta2_est, S.dist=S.dist, TT=seg2_tmp, T.burn=100))
      #Q.grid is where the failure comes, something in pl()
      q_emp[b] <- Q.grid(lambda.hat=seg1_tmp, TT=TT_tmp, S.dist=S.dist, y.full=yB,
                         theta1=theta1_est, theta2=theta2_est,
                         t.lag=t.lag, remedy=remedy, q.bound=round(0.2*TT_tmp)+1, ts=F,mu_zero=T)
    }
    t_compute3 <- proc.time() - t0
    
    ci <- quantile(q_emp, probs=c(0.025,0.975))
    ci_result <- rbind(ci_result, pelt_result[cp_index+1]+ci)
  }
}
num_cp
# save.image(paste0('Result_', data_name, '.RData'))
load('savResults.Rdata')
if (is.null(ci_result)) {
  savResults[[q]]<-list(NA,pelt_result)
  ci_result<-NA
  save(file = 'savResults.RData',list='savResults')
  save(list=c('ci_result','pelt_result'),file=paste(siz,q,'.RData',sep=''))
  return(list(NA,pelt_result))
} else {
  savResults[[q]]<-list(ci_result,pelt_result)
  save(file = 'savResults.RData',list='savResults')
  save(list=c('ci_result','pelt_result'),file=paste(siz,q,'.RData',sep=''))
  return(list(ci_result,pelt_result))
}
}
timEnd<-proc.time()
save(file=paste('FullCLMDL_',siz,'.RData',sep=''),
     list=c('savResults','timStart','timEnd'))
}