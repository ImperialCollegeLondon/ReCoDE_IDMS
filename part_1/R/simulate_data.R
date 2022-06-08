# Function to produce simulated reported incidence data -----------------------# 


# Input:

# - n_pop:    population size, assumed for Gauteng 
# - immunity: SARS-CoV-2 seroprevalance (%), user defined 
# - n_inf = 1: number initially infectious, assumed 1 
# - sigma = 1/5.1: progression rate, 1/ incubation period assumed to be 5.1 days  
# - gamma = 1/2.1: recovery rate, 1 / infectious period assumed to be 2.1 days 
# - rho  = reporting rate, used defined 
# - beta = transmission rate, user defined 
# - ts = time steps, user defined 

# Output: 

# - data frame of solutions to the derivatives of all compartments at each time step

simulate_data = function(
  n_pop = 15810388,               
  immunity , 
  n_inf ,
  sigma = 1/5.1,   
  gamma = 1/2.1,   
  rho ,
  beta , 
  ts
){
  
  # required package 
  library(deSolve)
  library(cowplot)

  # source model 
  source("models/model1_deSolve.R")
  
  # initial values 
  n_recov = n_pop * immunity 

  initial_state = c(S= n_pop - n_recov - n_inf  , E = 0 , I=n_inf ,Q=0, R=n_recov)
  params = c(S0 = n_pop - n_recov , beta=beta, sigma=sigma, gamma=gamma,rho=rho)
  
  model = ode(initial_state, ts, SEIQR, params)
  data  = round(as.data.frame(model))

  for(i in 2:ncol(data)) plot(data[,i], main = names(initial_state)[(i-1)] )
  
  
  return(data)
}

