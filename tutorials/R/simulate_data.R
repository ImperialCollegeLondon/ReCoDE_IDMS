# Function to produce simulated reported incidence data for single variant -----------------------# 


# Input:

# - n_pop:    population size, assumed for Gauteng (class = numeric)
# - immunity: SARS-CoV-2 seroprevalance (%), user defined (class = numeric, 0-1)
# - n_inf: number initially infectious (class = numeric)
# - sigma:  progression rate, 1/ latent period assumed to be 3.03 days (class = numeric)
# - gamma: recovery rate, 1 / infectious period assumed to be 4.17 days (class = numeric)
# - rho: reporting probability, used defined (class = numeric)
# - beta: transmission rate, user defined (class = numeric)
# - ts: vector of time steps, user defined (class = numeric)

# Output: 

# - data frame of solutions to the derivatives of all compartments at each time step

simulate_data_single_var = function(
  n_pop = 15810388,               
  immunity , 
  n_inf ,
  gamma = 1/4.17,                
  sigma = 1/3.03,   
  rho ,
  beta , 
  ts
){
  
  # required package 
  library(deSolve)


  # source model 
  source("models/model1_deSolve.R")
  
  if(n_pop <0 ) stop("n_pop cannot be negative")
  if(immunity < 0 | immunity > 1) stop("immunity should be between 0 and 1")
  if(rho < 0 | rho > 1) stop("Rho should be between 0 and 1")
  if(n_inf < 0) stop("n_inf cannot be negative")
  if(gamma < 0) stop("gamma cannot be negative")
  if(sigma < 0) stop("sigma cannot be negative")
  if(beta < 0) stop("beta cannot be negative")  
  if(length(ts) < 2) stop("must be at least two time steps")

  
  # initial values 
  n_recov = n_pop * immunity 

  initial_state = c(S= n_pop - n_recov - n_inf  , E = 0 , I=n_inf ,Q=0, R=n_recov)
  params = c(S0 = n_pop - n_recov , beta=beta, sigma=sigma, gamma=gamma,rho=rho)
  
  model = ode(initial_state, ts, SEIQR, params)
  data  = round(as.data.frame(model))

  for(i in 2:ncol(data)) plot(data[,i], main = names(initial_state)[(i-1)] )
  
  
  return(data)
}


# Function to produce simulated reported incidence data for 2 variants -----------------------# 


# Input:

# - n_pop:    population size, assumed for Gauteng (class = numeric)
# - immunity: SARS-CoV-2 seroprevalance (%), user defined (class = numeric, 0-1)
# - n_inf_D: number initially infectious, used defined (class = numeric)
# - n_inf_O: number initially infectious, used defined (class = numeric)
# - sigma_D:  progression rate, 1/ latent period, used defined  (class = numeric)
# - gamma_D: recovery rate, 1 / infectious period, used defined (class = numeric)
# - sigma_O:  progression rate, 1/ latent period, used defined  (class = numeric)
# - gamma_O: recovery rate, 1 / infectious period, used defined (class = numeric)
# - rho_D: reporting probability, used defined (class = numeric)
# - rho_O: reporting probability, used defined (class = numeric)
# - beta_D: transmission rate, user defined (class = numeric)
# - beta_O: transmission rate, user defined (class = numeric)
# - epsilon: immune waning rate, user defined (class = numeric)
# - omega: impact of interventions on transmission (%), user defined (class = numeric, 0-1)
# - nu: vaccination rate, assumed to be constant in time, user defined (class = numeric)
# - ts: vector of time steps, user defined (class = numeric)
# - time_int_start: time to start interventions, user defined  (class = numeric)
# - time_int_end: time to end interventions, user defined  (class = numeric)
# - time_seed_O: time to seed Omicron, user defined (class = numeric) 

# Output: 

# - data frame of solutions to the derivatives of all compartments at each time step

simulate_data_multi_var = function(
  n_pop = 15810388,               
  immunity , 
  n_inf_D ,
  n_inf_O ,
  gamma_D,                
  sigma_D,   
  gamma_O,               
  sigma_O, 
  rho_D ,
  rho_O,
  beta_D ,
  beta_O,
  epsilon, 
  omega,
  nu, 
  ts,
  time_int_start,
  time_int_end,
  time_seed_O
){
  
  # required package 
  library(deSolve)
  
  
  # source model 
  source("models/model2_deSolve.R")
  
  
  
  if(n_pop <0 ) stop("n_pop cannot be negative")
  if(immunity < 0 | immunity > 1) stop("immunity should be between 0 and 1")
  if(rho_D < 0 | rho_D > 1 | rho_O < 0 | rho_O > 1) stop("rho should be between 0 and 1")
  if(omega < 0 | omega > 1) stop("omega should be between 0 and 1")
  if(n_inf_D < 0 | n_inf_O < 0) stop("n_inf cannot be negative")
  if(gamma_D < 0 | gamma_O < 0) stop("gamma cnanot be negative")
  if(sigma_D < 0 | sigma_O < 0) stop("sigma cannot be negative")
  if(beta_D < 0 | beta_O < 0) stop("beta cannot be negative") 
  if(epsilon <0 ) stop("epsilon cannot be negative")
  if(nu <0 ) stop("nu cannot be negative")
  if(length(ts) < 2) stop("must be at least two time steps")
  if(time_int_start < ts[1] | time_int_start > tail(ts,1) | time_int_end < ts[1] | time_int_end > tail(ts,1))
    stop("interventions must be during modelling time period")
  if(time_int_start > time_int_end) stop ("interventions must start before they can end")
  if( time_seed_O < ts[1] | time_seed_O > tail(ts,1)  ) stop("omicron seed must be during modelling time period")
  

  
  # initial values 
  n_recov = n_pop * immunity 
  
  
  initial_state = c(
    S  = n_pop - n_recov - n_inf_D,
    ED = 0 ,
    ID = n_inf_D ,
    EO = 0,
    IO = 0,
    QD = 0,
    QO = 0,
    R  = n_recov,
    SO = 0
  )
  
  
  
  params = c(
    betaD = beta_D,
    betaO = beta_O,
    sigmaD = sigma_D,
    gammaD = gamma_D,
    sigmaO = sigma_O,
    gammaO = gamma_O,
    rhoD = rho_D,
    rhoO = rho_O, 
    epsilon = epsilon, 
    omega = omega,
    time_switch = time_int_start,
    time_switch_off =time_int_end,
    nu = nu,
    N = n_pop - n_recov 
  )
  
  
  # We use deSolve events to trigger the seeding of Omicron during the model simulation 
  seed_omicron_event =  data.frame(var = c("IO", "S"),
                                   time = time_seed_O,
                                   value = c(n_inf_O,-n_inf_O),
                                   method = c("add" ,"add"))
  
  model = data.frame(ode(initial_state, ts, SEIQR2, params, events = list(data = seed_omicron_event)))
  

  data  = round(as.data.frame(model))
  
  for(i in 2:ncol(data)) plot(data[,i], main = names(initial_state)[(i-1)] )
  
  
  
  return(data)
}


