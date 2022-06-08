# Function to calculate simulated reported incidence --------------------------# 


# Input:

# - ODE_data: solution to the derivatives of all compartments at each time step,
#   user defined (class = data.frame)
# - all_dates: list of all dates, user defined (class = dates)
# - date_seed: date to seed variant, user defined (class = date)
# - rho: reporting probability, user defined (class = numeric, 0-1)
# - sigma: recovery rate, assumed to be 1/5.1 (optional, class = numeric)

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

  data = data.frame(date = all_dates,
                        rep_inc = round(rho *ODE_data$E * sigma))
  
  data = filter(data, date >= date_seed)
  
  plot = ggplot(data, aes(x= date , y = rep_inc)) +
    geom_point() +
    ylab("reported incidence") +
    xlab("date")
  

  
  return(list(data, plot))
}
