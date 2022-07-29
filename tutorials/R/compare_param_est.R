# Function to compare parameter estimates between models or methods of --------#
# solving ODEs.  --------------------------------------------------------------# 


# Input:

# - parameter_names: names of parameters to compare, user defined (class = character).
# - true_param_values: vector of true parameter values, user defined (class = numeric)
# - param_values1: mean and 95% CrI of parameter values, user defined estimated by model 1
#   first column = mean, second column = 2.5% CrI, third column = 97.5% column (class = matrix or data.frame). 
# - param_values2: mean and 95% CrI of parameter values, user defined estimated by model 2 
#   first column = mean, second column = 2.5% CrI, third column = 97.5% column(class = matrix or data.frame). 
# - model_names: if comparing names of the  model(s) / method(s) being compared, assumed to be EU and RK (class = character)

# Output: 

# - list of ggplots

compare_param_est = function(
  parameter_names,
  true_param_values,
  param_values1,
  param_values2 = NULL,  
  model_names = c("EU", "RK")
){
  
  if(length(model_names)==2 & is.null(param_values2)) stop ("if comparing 2 models, 2 sets of estimate parameter values required")
  if(! is.data.frame(param_values1) &   ! is.matrix(param_values1)) stop ("param_values1 must be a date.frame or matrix" )
  if(! is.null(param_values2) & is.data.frame(param_values1) &   ! is.matrix(param_values1)) stop ("param_values2 must be a date.frame or matrix" )

    
  # required package 
  library(ggplot2)
  library(dplyr)
  
  true_param = data.frame(names = parameter_names,
                          mean = true_param_values)
  
  param1 = data.frame(names = parameter_names, 
                      mean  =  param_values1[,1],
                      lower =  param_values1[,2],
                      upper =  param_values1[,3])
  
  if(length(model_names)==2 & ! is.null(param_values2) ){
    
    param2 = data.frame(names = parameter_names, 
                        mean  = param_values2[,1],
                        lower = param_values2[,2],
                        upper = param_values2[,3])
  

  param_df = true_param  %>% 
    bind_rows(param1) %>% 
    bind_rows(param2) %>%  
    mutate(model = c(rep("true", length(parameter_names)), 
                     rep(model_names[1], length(parameter_names)), 
                     rep(model_names[2],length(parameter_names)))) 
  } else if (length(model_names) == 1) {
    
    param_df =  true_param %>% 
      bind_rows(param1) %>% 
      mutate(model = c(rep("true", length(parameter_names)), 
                       rep(model_names[1], length(parameter_names)))) 
  } else {
    stop("Can only compare up to 2 models to true values")
  }
  
  param_plots = list()
  for(i in 1:length(parameter_names)){
    param_plots[[i]] =  param_df %>% 
      filter(names == parameter_names[[i]]) %>% 
      ggplot(aes(x=names, y = mean))+
      geom_point(aes(color = model), position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = lower, ymax = upper, color = model),
                    position = position_dodge(width = 0.5), width =  0.4, size = 1) +
      ylab("parameter estimate") 
  }
  
  
  
  return(param_plots)
}

