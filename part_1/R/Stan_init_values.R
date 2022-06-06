# Function to generate a plausible range of parameter values ------------------# 
# for each chain to start at during Rstan model fitting. ----------------------#  

# Input: 

# - seed for reproducability 

# Output: 

# - starting values for our 3 model parameters, beta, I0, rho 

ini_1 = function(seed = 1){
  set.seed(seed)
  list(beta=runif(1,1.7,5.2),
       I0 = runif(1,1,200),
       rho = runif(1,0,0.5) )}
