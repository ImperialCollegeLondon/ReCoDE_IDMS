# Function to plot the results of a singel variant fitted Stan model ----------#
# against the data to which it was fit  ---------------------------------------# 


# Input:

# - stan_fit: model fit results, user defined (class = stanfit).
# - variable_model: name of the model variable to plot, user defined (class = character)
# - variable_data: name of the data variable to plot, user defined (class = character). 
# - data: observed or simulated data which the model was fit to. One column named "date" 
# - should be of class dates,  user defined (class = data.frame). 
# - all_dates: vector of dates over which the model was run (class = date)

# Output: 

# - ggplot of the model fit 

plot_model_fit = function(
 stan_fit,
 variable_model,
 variable_data,
 data, 
 all_dates
){
  
  
  # required package 
  library(ggplot2)
  library(dplyr)
  
  # extract posterior estimates from stan fit 
  stan_fit_ext = rstan::extract(stan_fit) 
  
  stan_fit_df = stan_fit_ext[[variable_model]] %>% as.data.frame.table() 
  
  # group by data and calculate the mean and 95% CrI 
  stan_fit_df_sum = stan_fit_df %>% 
    rename(ni = iterations, time = Var2, value = Freq) %>%
    mutate(
      ni = as.numeric(ni),
      time = as.numeric(time)) %>% 
    group_by(time) %>% 
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) %>%  
    mutate(date= all_dates)
  
  
  ## create a name from the character string describing the variable to plot 
  variable_data = sym(variable_data)
  
  # left join observed data and plot using ggplot 
  plot = stan_fit_df_sum %>% 
    left_join(data, by = "date") %>%  
    ggplot(aes(x = date , y = !!variable_data)) +
    geom_point()+
    geom_line(aes(y=mean)) +
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha =0.2) +
    ylab("Reported Incidence") + xlab("Date")
  
  
  return(plot)
}


# Function to plot the results of a fitted Stan model with 2 varaints ---------#
# against the data to which it was fit  ---------------------------------------# 


# Input:

# - stan_fit: model fit results, user defined (class = stanfit).
# - variable_model: name of the model variable to plot, user defined (class = character)
# - variable_data: name of one of the data variable to plot, should correspond to first named variant (see variant_names) user defined (class = character). 
# - data: observed or simulated data which the model was fit to, user defined (class = data.frame).
# - variant_names: vector of names of variants being plotted, assumed to be Delta and Omicron (class = character)

# Output: 

# - ggplot of the model fit 

plot_model_fit_multi_var = function(
  stan_fit,
  variable_model,
  variable_data,
  data,
  variant_names = c("Delta","Omicron")
){
  
  
  # required package 
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # extract posterior estimates from stan fit 
  stan_fit_ext = rstan::extract(stan_fit) 
  
  stan_fit_df = stan_fit_ext[[variable_model]] %>% as.data.frame.table() 
  
  # group by data and calculate the mean and 95% CrI 
  stan_fit_df_sum = stan_fit_df %>% 
    rename(ni = iterations, time = Var2, variant = Var3, value = Freq) %>%
    mutate(
      ni = as.numeric(ni),
      time = as.numeric(time),
      variant = ifelse(variant == "A", variant_names[1], variant_names[2])) %>% 
    group_by(time,variant) %>% 
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) 
  
  # here we make a variable from the row names of our data, 
  # which account for the discarded unobserved data before fitting

  data2 =  data %>% 
    mutate(time = (1:dim(data)[1])) %>%  
    pivot_longer(cols = c(rep_inc_D_noise, rep_inc_O_noise)) %>%
    mutate(variant = ifelse(name == variable_data, variant_names[1], variant_names[2])) %>% 
    select(date, time, value, variant)
  

  # left join observed data and plot using ggplot 
  plot = stan_fit_df_sum %>% 
    left_join(data2) %>%  
    ggplot(aes(x = time , y = value)) +
    geom_point(aes(color=variant))+
    geom_line(aes(y=mean, color = variant)) +
    geom_ribbon(aes(ymin=lower,ymax=upper, fill = variant), alpha =0.2) 
  
  return(plot)
}

