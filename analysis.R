rm(list=ls())
library(trend)

data_name <- 'data_n450' # one change-point
load(paste0(data_name, ".RData"))

y <- data$y # T times S spatio-temporal data
S.dist <- data$S.dist # distance matrix of the S spatial locations
S <- ncol(y) # spatial dimension
TT <- nrow(y) # time dimension

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

# Summarize results
summary <- data.frame(
  Location = 1:ncol(y),
  Change_Point = change_points,
  P_Value = p_values
)

colnames(summary,)
print(summary)

