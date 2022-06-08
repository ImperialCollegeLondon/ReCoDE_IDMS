# Function to run diagnostics on a stan fit -----------------------------------# 


# Input:

# - stan_fit: an object of S4 class which contains the model fit results.
# - pars: names of parameters we want to check 

# Output: 

# - Diagnostic plots
# - data frame of summary statistics 

diagnose_stan_fit = function(
stan_fit,
pars){

  # required package 
  library(bayesplot)
  library(cowplot)


  stan_fit_post= as.array(stan_fit)

  # markov chain trace plots   
  markov_trace = mcmc_trace(stan_fit_post, pars=c("lp__", pars))
  
  # bivariate marginal posterior distributions
  pairs(stan_fit, pars = pars, cex.labels=1.5, font.labels=9, condition = "accept_stat__")  
  
  # univariate marginal posterior distributions 
  uni = mcmc_dens_overlay(stan_fit, pars=pars)
  
  # summary of parameter values, effective sample size and Rhat 
  param_sum = summary(stan_fit, pars = pars)$summary

  # output files 
  list_out= list( markov_trace,uni, param_sum)
  names(list_out) = c("markov chain trace plots",
                      "univariate marginal posterior distributions", "summary statistics of parameters")
 return(list_out)

}
