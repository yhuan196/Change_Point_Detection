# Project: Change_Point_Detection

This folder contains 3 files, utilities_generate.R, utilities_estimate.R, utilities_main.R.


-------------------------------------------------------------------------------------------------------------
## utilities_generate.R
This file contains the R code to generate a spatio-temporal process data on four autoregression spatial model.

 1. `Matern()`: is the funtion that produce Matern covariance.
 2. `sim.()`: funtion simulates one spatio-temporal segment.
 3. `gen_one_change_point_data()`: generates one-change-point spatio-temporal data.
-------------------------------------------------------------------------------------------------------------
 

-------------------------------------------------------------------------------------------------------------
## utilities_estimate.R
This file contains all the functions to perform CLMDL analysis.

  1. `s.dist.pair()` and `s.dist.pair.0()`: Auxiliary function to compute pairs within s.lag distance.
  3. `D.cal()`: Input y and distance matrix find the D.pair (The total number of pairs of observations within 
  the specified distance.), use.obs (the number of times each observation is used in pair calculations) and 
  a remedy matrix for the edge effect.
  4. `pl()`: Calculate -logPL and D.pair with given distance matrix, lags and remedy matrix.
  5. `cost()` Cost function for a segment.
  6. `pelt()`: Function that outputs the change-point position by considering the whole dataset.
  7. `total_cost`: Given change point location, calculate CLMDL (cost function).
  8. `Q.grid`: function to calculate CI.
-------------------------------------------------------------------------------------------------------------


-------------------------------------------------------------------------------------------------------------
## utilities_main.R
This is the file to implement Pettitt's Test and CLMDL algorithm on the four-parameter autoregressive spatial 
model. 
  
  1. `run_simulation()` is the simulation on the four-parameter autoregressive spatial model using Pettitt's 
  Test. The output is summary statistics including estimated Lambda, standard deviation of Lambda, coverage 
  proportion, 95% confidence interval, and detection rate.
  
  2. The final result of CLMDL algorithm `savResults` is from the last loop. In the loop, we use the the `pelt()` 
  function detects change points, then estimate the Lambda for segments, and employs a bootstrap to construct 
  confidence intervals for change point locations.
-------------------------------------------------------------------------------------------------------------
  