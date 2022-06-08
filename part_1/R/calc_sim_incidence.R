# Function to calculate simulated reported incidence --------------------------# 


# Input:

# - ODE_data: solution to the derivatives of all compartments at each time step (class = data.frame)
# - all_dates: list of all dates (class = dates)
# - date_seed: date to seed variant (class = character)
# - rho: reporting rate (class = numeric)
# - sigma: recovery rate (optional, class = numeric)

# Output: 

# List of: 
# - date frame of simulated reported incidence over time 
# - ggplot figure of simulated reported incidence over time

calc_sim_incidence = function(
  ODE_data,
  all_dates,
  date_seed,
  rho,
  sigma = 1/5.1
){
  
  # required package 
  library(tidyverse)

  data = data.frame(time = all_dates,
                        rep_inc = round(rho *ODE_data$E * sigma))[date_seed:length(all_dates), ]
  
  plot = ggplot(data, aes(x= time , y = rep_inc)) +
    geom_point() +
    ylab("reported incidence") +
    xlab("date")
  

  
  return(list(data, plot))
}
