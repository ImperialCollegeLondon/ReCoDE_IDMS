
# Function to solve 2 variant SEIQR model -----------------------------------------------#  

# Input: 

# - time: time steps to fit model, used defined (class = numeric) 
# - current_state: initial values of all compartments, used defined (class = numeric)
# - params: rate parameters and time points at which to introduce interventions, used defined (class = numeric)

# Output: 

# - derivatives of all compartments at each time step 


SEIQR2 = function(times, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S + ID + ED + EO + IO + R + QD + QO + SO
    
    ##  model interventions

    
    if(times < time_switch){ 
      foiD =  betaD*(ID/N)
      foiO =  betaO*(IO/N)
    } else if (times > time_switch_off){ 
      foiD =  betaD*(ID/N)
      foiO =  betaO*(IO/N)
    }  else {
      foiD = (1-omega) * betaD*(ID/N)
      foiO = (1-omega) * betaO*(IO/N)
    }
  
    
    
    dS = - (foiD+foiO)* S - nu * S 
    dED = foiD * S - sigmaD * ED
    dEO = foiO * (S+SO) - sigmaO * EO
    
    dID = (1-rhoD) * sigmaD * ED - gammaD * ID 
    dIO = (1-rhoO) * sigmaO * EO - gammaO * IO 
    
    dQD =  rhoD* sigmaD * ED  - gammaD * QD
    dQO =  rhoO * sigmaO * EO - gammaO * QO
      
    dR = gammaD * (ID+ QD) + gammaO * (IO + QD) - epsilon * R + nu * S
    
    dSO = epsilon * R - foiO * SO 
    
    
    return(list(c(dS,dED, dID,dEO, dIO,dQD,dQO, dR, dSO)))
  })
}
