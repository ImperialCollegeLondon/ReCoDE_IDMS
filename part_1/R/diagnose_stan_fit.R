# Function to run diagnostics on a stan fit -----------------------------------# 


# Input:

# - stan_fit: an object of S4 class which contains the model fit results.
# - pars: names of parameters we want to check 

# Output: 

# - Number of divergent transitions 
# - Diagnostic plots
# - data.frame of summary statistics 

diagnose_stan_fit = function(
stan_fit){

  # required package 
  library(bayesplot)
  library(cowplot)
  
  # check divergent transitions: 
  divergent = check_divergences(stan_fit)
  

  stan_fit_post= as.array(stan_fit)

  # markov chain trace plots   
  markov_trace = mcmc_trace(stan_fit_post, pars=c("lp__", pars))
  
  # bivariate marginal posterior distributions
  bivar = pairs(m1_fit_euler, pars = pars, cex.labels=1.5, font.labels=9, condition = "accept_stat__")  
  
  # univariate marginal posterior distributions 
  uni = mcmc_dens_overlay(m1_euler_post, pars=pars)
  
  plot_out = plot_grid(markov_trace, bivar, uni, ncols = 3)

  # summary of parameter values, effective sample size and Rhat 
  param_sum = summary(m1_fit_euler, pars = pars)$summary

 return(list(divergent, plot_out, param_sum))

}
