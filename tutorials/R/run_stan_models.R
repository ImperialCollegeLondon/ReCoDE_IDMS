# Function to fit a stan model ------------------------------------------------# 


# Input:

# - list_data: data for model fitting, user defined (class = list)
# - model: model compiled from C++ code, user defined (class = stanmodel)
# - n_chains: number of Markov chains, assumed to be 3 (class = numeric, positive integer)
# - n_iter: number of iterations for each chain, assumed to be 2000 
#   (class = numeric, positive integer)
# - n_warmup: number of warmup iterations per chain, assumed to be 1000 
#   (class = numeric, positive integer, must be less than n_iter)
# - seed_values: vector of values for setting seed when choosing initial values, 
#   assumed to be 1,2,3 (class = numeric, positive integer, 
#   vector length must equal n_chains)
# - n_var: number of variants, assumed to be 1 (class = numeric)
# - model_no: model number to fit (1= 1 variant, 2 = 2 variant), determines which parameters to draw initial values for (class = numeric)

# Output: 

# - an object of S4 class which contains the model fit results.

run_stan_models = function(
  list_data, 
  model,
  n_chains = 3,
  n_iter = 2000,
  n_warmup = 1000, 
  seed_values = c(1,2,3),
  n_var =1,
  model_no=1){

  # required packages
  library(rstan)
  rstan_options(auto_write = TRUE)           
  options(mc.cores = parallel::detectCores())
  
  
  if(length(seed_values) != n_chains){
    stop("seed values do not equal number of chain")
  }

  
  # source function to draw initial values 
source("R/draw_init_values.R")
  
  list_of_inits = list()

  if(model_no == 1){  
    for(i in 1:n_chains)  {
      list_of_inits[[i]] =  draw_init_values(seed = seed_values[i],
                                             n_var = n_var
      )}} else {
        for(i in 1:n_chains)  {
          list_of_inits[[i]] =  draw_init_values_2(seed = seed_values[i],
                                                 n_var = n_var
          )}
      }
  
  
  
  # draw samples from Stan model 
  
time.start = Sys.time() 
  stan_fit = sampling(
    model,
    data = list_data,
    chains=n_chains,
    warmup=n_warmup,
    iter=n_iter,
    init = list_of_inits,
    control = list(
      max_treedepth = 15,
      adapt_delta= 0.87
    )
  )
 time.end = Sys.time()
 
 print(time.end - time.start)
 
 return(stan_fit)
}
