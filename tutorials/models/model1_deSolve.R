
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
