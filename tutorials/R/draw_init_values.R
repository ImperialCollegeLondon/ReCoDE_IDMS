# Function to generate a different starting values for each Markov chain ------# 
# Single variant model --------------------------------------------------------#

# Input:  

# - seed: assumed to be 1 (class = numeric)
# - n_var: number of variants (class = numeric)


# Output: 

# - starting values for our model parameters, beta, I0, rho, k 

draw_init_values = function(seed = 1,
                            n_var = 1){
  set.seed(seed)
 list(
       beta=runif(n_var,1,4),
       rho = runif(n_var,0,.5),
       k = runif(1,0,1))
  }


draw_init_values_2 = function(seed = 1,
                              n_var = 1){
  set.seed(seed)
  list(
    k = runif(1,0,5),
    epsilon = runif(1, 0,0.011)
    )
}


