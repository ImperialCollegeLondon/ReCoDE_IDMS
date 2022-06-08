# Function to plot the results of a fitted Stan model -------------------------#
# against the data to which it was fit  ---------------------------------------# 


# Input:

# - stan_fit: an object of S4 class which contains the model fit results.
# - variable_model: name of the model variable to plot
# - variable_data: name of the data variable to plot. 
# - data: observed or simulated data which the model was fit to. 

# Output: 

# - ggplot of the model fit 

plot_model_fit = function(
 stan_fit,
 variable_model,
 variable_data
 data 
){
  
  # required package 
  library(tidyverse)
  
  # extract posterior estimates from stan fit 
  stan_fit_ext = rstan::extract(stan_fit) 
  
  stan_fit_df = stan_fit_ext$variable %>% as.data.frame.table() 
  
  # group by data and calculate the mean and 95% CrI 
  post_euler_df_sum = post_euler_df %>% 
    rename(ni = iterations, time = Var2, value = Freq) %>%
    dplyr::mutate(
      ni = as.numeric(ni),
      time = as.numeric(time)) %>% 
    group_by(time) %>% 
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) 
  
  # left join observed data and plot using ggplot 
  
 plot = post_euler_df_sum %>% 
  left_join(sim_data) %>%  
    ggplot(aes(x = time , y = variable_data)) +
    geom_point()+
    geom_line(aes(y=mean)) +
    geom_ribbon(aes(ymin=lower,ymax=upper)) 
  
  return(plot)
}

