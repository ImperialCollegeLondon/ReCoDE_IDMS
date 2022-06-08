# Function to produce simulated reported incidence data -----------------------# 


# Input:

# - n_pop:    population size, assumed for Gauteng (class = numeric)
# - immunity: SARS-CoV-2 seroprevalance (%), user defined (class = numeric, 0-1)
# - n_inf: number initially infectious, assumed 1 (class = numeric)
# - sigma:  progression rate, 1/ incubation period assumed to be 5.1 days (class = numeric)
# - gamma: recovery rate, 1 / infectious period assumed to be 2.1 days (class = numeric)
# - rho: reporting probability, used defined (class = numeric)
# - beta: transmission rate, user defined (class = numeric)
# - ts: vector of time steps, user defined (class = numeric)

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

