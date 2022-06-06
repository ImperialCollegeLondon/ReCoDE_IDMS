# Function to produce simulated reported incidence data -----------------------# 


# No input:

# - n_pop:    population size, assumed for Gauteng 
# - immunity: SARS-CoV-2 seroprevalance (%), user defined 
# - n_inf = 1: number initially infectious, assumed 1 
# - sigma = 1/5.1: progression rate, 1/ incubation period assumed to be 5.1 days  
# - gamma = 1/2.1: recovery rate, 1 / infectious period assumed to be 2.1 days 
# - rho  = reporting rate, used defined 
# - beta = transmission rate, user defined 
# - ts = time steps, user defined 

# Output: 

# - starting values for our 3 model parameters, beta, I0, rho 

simulate_data = function(
  n_pop = 15810388                
  immunity , 
  n_inf = 1,
  sigma = 1/5.1   
  gamma = 1/2.1   
  rho ,
  beta , 
  ts
){
  
  # required package 
  library(deSolve)

  # source model 
  source("models/model1_deSolve")
  
  # initial values 
  n_recov = n_pop * immunity 
  S0 = n_pop - n_recov 
  
  
  initial_state = c(S= S0 - n_inf  , E = 0 , I=n_inf ,Q=0, R=n_recov)
  params = c(beta, sigma, gamma,rho,S0)
  
  model = ode(initial_state, ts, SEIQR, params)
  out.df  = round(as.data.frame(model))
  
  return(out.df)
}

