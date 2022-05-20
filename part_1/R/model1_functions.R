##########################################
# Function used to simulate data 
##########################################

# Function to solve SEIQR model

# Input: 

# - time: time steps to fit model
# - current_state: initial values of all compartments
# - params: rate parameters 

# Output: 

# - derivatives of all compartments at each time step 


SEIQR = function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N = S+I+E+Q+R
    
    dS = -(beta*S*I)/N
    
    dE = (beta*S*I)/N - sigma * E
    
    dI = (1-rho) * sigma * E - gamma*I 
    
    dQ = rho * sigma * E - gamma * I 
    
    dR = gamma*(I+Q)
    
    
    return(list(c(dS,dE, dI,dQ, dR)))
  })
}
