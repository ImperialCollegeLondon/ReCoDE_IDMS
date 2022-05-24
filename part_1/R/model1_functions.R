################################################################################
# This script contains function used in part 1                                 #
################################################################################



# Function to solve SEIQR model -----------------------------------------------#  

# Input: 

# - time: time steps to fit model
# - current_state: initial values of all compartments
# - params: rate parameters 

# Output: 

# - derivatives of all compartments at each time step 



SEIQR = function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{

    dS <- -(beta*S*I)/S0
    dE <- (beta*S*I)/S0 - sigma*E
    dI <- (1-rho)*sigma*E - gamma*I 
    dQ <- rho * sigma * E - gamma * Q
    dR <- gamma*(I+Q)
    
    return(list(c(dS, dE, dI,dQ, dR)))
  })
}


# Function to generate a plausible range of parameter values ------------------# 
# for each chain to start at during Rstan model fitting. ----------------------#  

# No input, the function just samples randomly from a range of values

# Output: 

# - starting values for our 3 model parameters, beta, I0, rho 

ini_1 = function(seed = 1){
  set.seed(seed)
  list(beta=runif(1,1.7,5.2),
       I0 = runif(1,1,200),
       rho = runif(1,0,0.5) )}
