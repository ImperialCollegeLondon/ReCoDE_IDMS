
# Function to tidy raw reported incidence data -----------------------# 

# Input:

# - raw_data: cumualtive reported incidence over time, user defined (class = data.frame,
# - minimum columns should be date and incidence varaible)
# - date_seed: date to seed variant, user defined (class = date)
# - last_date: last date of modelling period, user defined (class = date)
# - var: name of incidence variable, assumed to be GP (class = character)


# Output: 

# - data frame of daily reported incidence and cumulative incidence over the 
#   specified time period. 

tidy_obs_inc = function(
  raw_data,
  date_seed,
  end_date,
  var = "GP"
){
  
  library(tidyverse)
  
   var = sym(var)
  
   data_out = raw_data %>%
    select(date, var) %>%
    rename(cum_rep_inc = var) %>%
    # extract daily incidence from cumulative 
    mutate(rep_inc =  c(cum_rep_inc[1], diff(cum_rep_inc))) %>% 
    mutate(date = as.Date.character(date, format = "%d-%m-%Y")) %>%
    filter(date >= date_seed & date <= end_date)
   
   
   return(data_out)
}

