# Function to generate a different starting values for each Markov chain ------# 


# Input:  

# - seed, assumed to be 1 (class = numeric)

# Output: 

# - starting values for our 3 model parameters, beta, I0, rho 

draw_init_values = function(seed = 1){
  set.seed(seed)
  list(beta=runif(1,1,12),
       I0 = runif(1,1,200),
       rho = runif(1,0,1) )}
