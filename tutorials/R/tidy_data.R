# Function to calculate a moving average -----------------------------#


# Input:

# - data: incidence over time, user defined (class = numeric)
# - time_int: time period over which to average (class = numeric) 

# Output:

# - vector of moving average incidence 

moving_av = function(data, time_int){stats::filter(data, rep(1 /time_int,time_int), sides = 2)}

# Function to tidy raw reported incidence data -----------------------# 

# Input:

# - raw_data: cumulative reported incidence over time, user defined (class = data.frame)
# - minimum columns should be date and incidence variable)
# - date_seed: date to seed variant, user defined (class = date)
# - last_date: last date of modelling period, user defined (class = date)
# - var: name of incidence variable, assumed to be GP (Gauteng province) (class = character)
# - time_int: time period over which to take a moving average, assumed to be 7 (class = numeric) 



# Output: 

# - data frame of the date, the cumulative incidence, the daily reported incidence 
#   and moving average of the daily reported incidence, over the specified time period. 

tidy_obs_inc = function(
  raw_data,
  date_seed,
  end_date,
  var = "GP",
  time_int=7
){
  
   # load package
  
  library(dplyr)
  
  
  all_dates = data.frame(
    date = seq.Date(from = date_seed, to = end_date ,  by = "days")) # model times
  
  
   var = sym(var)
  
   data_tidy = raw_data %>%
    select(date, var) %>%
    rename(cum_rep_inc = var) %>%
    # extract daily incidence from cumulative 
    mutate(rep_inc =  c(cum_rep_inc[1], diff(cum_rep_inc)),
           rep_inc_ma = round(moving_av(rep_inc, time_int=time_int)),
           date = as.Date.character(date, format = "%d-%m-%Y")) %>%
    filter(date >= date_seed & date <= end_date)
   
   
   data_out = all_dates %>%   # in case of missing dates 
     left_join(data_tidy) %>%  
     mutate(rep_inc = ifelse(is.na(rep_inc), 0 , rep_inc),
            rep_inc_ma = ifelse(is.na(rep_inc_ma), 0 , rep_inc_ma))  # remove NA values
   
   
   return(data_out)
}


# Function to tidy raw vaccination data -----------------------# 

# Input:

# - raw_data: cumulative number of vaccinated over time, user defined (class = data.frame)
# - date_seed: date to seed variant, user defined (class = date)
# - last_date: last date of modelling period, user defined (class = date)
# - var: name of vaccination variable, assumed to be GP (Gauteng province) (class = character)


# Output: 

# - data frame of the date, the cumulative number vaccinated, the daily number vaccinated
#   and the per capita rate of vaccination, over the specified time period. 


tidy_vacc = function(
  raw_data,
  date_seed,
  end_date,
  var = "GP",
  n_pop 
){
  
  # load package
  
  library(dplyr)
  

  all_dates = data.frame(
    date = seq.Date(from = date_seed, to = end_date ,  by = "days")) # model times
  
  var = sym(var)
  
  data_tidy = raw_data %>% 
    select(date, var) %>%
    rename(cum_vac = var) %>%
    # extract daily incidence from cumulative 
    mutate(vac =  c(cum_vac[1], diff(cum_vac)),
           date = as.Date.character(date, format = "%Y-%m-%d"),
           vac_rate = vac / n_pop) %>%                                  # calc per capita rate of vaccination  
    filter(date >= date_seed & date <= end_date)
  
    data_out = all_dates %>%   # in case of missing dates 
      left_join(data_tidy) %>%  
      mutate(vac_rate = ifelse(is.na(vac_rate), 0 , vac_rate))  # remove NA values
    
    
  return(data_out)
}

