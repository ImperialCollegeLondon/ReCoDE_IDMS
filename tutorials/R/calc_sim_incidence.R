# Function to calculate simulated reported incidence for single variant --------------------------# 


# Input:

# - ODE_data: solution to the derivatives of all compartments at each time step,
#   user defined (class = data.frame)
# - all_dates: list of all dates, user defined (class = dates)
# - date_seed: date to seed variant, user defined (class = date)
# - rho: reporting probability, user defined (class = numeric, 0-1)
# - sigma: latent rate, assumed to be 1/3.03 (optional, class = numeric)

# Output: 

# List of: 
# - date frame of simulated reported incidence over time 
# - ggplot figure of simulated reported incidence over time

calc_sim_incidence_single_var = function(
  ODE_data,
  all_dates,
  date_fit,
  rho,
  sigma = 1/3.03
){
  
  # required package 
  library(ggplot2)
  library(dplyr)
  
  # calculate reported incidence
  data = data.frame(date = all_dates,
                        rep_inc = round(rho *ODE_data$E * sigma))

  
  # extract dates to fit 
  data = filter(data, date >= date_fit)

  # add noise to data
  set.seed(2)
  for (i in 1:nrow(data)){
    data$rep_inc_noise[i] = rnbinom(1, mu = data$rep_inc[i], size = 10)
   
  }
  
  plot = ggplot(data, aes(x= date , y = rep_inc)) +
    geom_line(aes(color = "Reported incidence")) +
    geom_line(aes(y = rep_inc_noise, color = "Reported incidence + noise"))+
    ylab("Reported incidence") +
    xlab("Date") + 
    scale_colour_manual(name = "",
                        values = c("Reported incidence" = "red",
                                   "Reported incidence + noise" = "black"))
  

  
  return(list(data, plot))
}



# Function to calculate simulated reported incidence for 2 variants --------------------------# 


# Input:

# - ODE_data: solution to the derivatives of all compartments at each time step, user defined (class = data.frame)
# - all_dates: list of all dates, user defined (class = dates)
# - date_fit_D: date to fit variant, user defined (class = date)
# - date_fit_O: date to fit variant, user defined (class = date)
# - rho_D: variant-specific reporting probability, user defined (class = numeric, 0-1)
# - rho_O: variant-specific reporting probability, user defined (class = numeric, 0-1)
# - sigma_D: latent rate, user defined (class = numeric)
# - sigma_O: latent rate, user defined (class = numeric)

# Output: 

# List of: 
# - data frame of simulated reported incidence over time 
# - ggplot figure of simulated reported incidence over time

calc_sim_incidence_multi_var = function(
  ODE_data,
  all_dates,
  date_fit_D,
  date_fit_O,
  rho_D,
  rho_O,
  sigma_D,
  sigma_O
){
  
  # required package 
  library(ggplot2)
  library(dplyr)
  
  # calculate incidence 
  data = data.frame(date = all_dates,
                    rep_inc_D = round(rho_D *ODE_data$ED * sigma_D),
                    rep_inc_O = round(rho_O *ODE_data$EO * sigma_O))
  
  # add noise to data 
  set.seed(14)
  for (i in 1:nrow(data)){
    data$rep_inc_D_noise[i] = rnbinom(1, mu = data$rep_inc_D[i], size = 10)
    data$rep_inc_O_noise[i] = rnbinom(1, mu = data$rep_inc_O[i], size = 10)
  }

  # extract variant specific dates to fit 
  data2= data %>%  
    mutate(rep_inc_D = ifelse(date >= date_fit_D, rep_inc_D, NA ),
           rep_inc_O = ifelse(date >= date_fit_O, rep_inc_O, NA ),
           rep_inc_D_noise = ifelse(date >= date_fit_D, rep_inc_D_noise, NA ),
           rep_inc_O_noise = ifelse(date >= date_fit_O, rep_inc_O_noise, NA ))
  
  plot = ggplot(data2, aes(x= date , y = rep_inc_D_noise)) +
    geom_line(aes(color = "Delta")) +
    geom_line(aes(y= rep_inc_O_noise, color = "Omicron")) + 
    ylab("reported incidence") +
    xlab("date") + 
    scale_color_manual(name='Variant',
                         values=c('Delta'="#F8766D", 'Omicron'="#00BFC4" ))
  
  
  
  return(list(data2, plot))
}

