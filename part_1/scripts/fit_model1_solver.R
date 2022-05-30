################################################################################
# This code fits model 1 tp the simulated reported incidence data for the      #
# emergence of Omicron in Gauteng, South Africa. Model fitting is done using   #
# Rstan and its in built ODE solver.                                           #
################################################################################



# set up ----------------------------------------------------------------------#  

rm(list = ls())
file_path = "part_1/"


# load packages and source user defined functions -----------------------------#       

library(rstan)
library(bayesplot)
library(loo)
library(tidyverse)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

source(paste0(file_path,"R/model1_functions.R"))

# read in simulated data to fit the model to ----------------------------------#  

sim_data = read.csv(paste0(file_path, "data/sim_data.csv"))

# define global variables -----------------------------------------------------#  

# start and end data of model fitting 
start_date =  as.Date.character("01-09-2021", format = "%d-%m-%Y") 
end_date = as.Date.character("23-02-2022", format = "%d-%m-%Y")

# sequence of dates of model fitting 
all_dates = seq.Date(from = start_date, to = end_date ,  by = "days")

# date to seed omicron 
seed_omicron =  which(all_dates == as.Date.character("01-10-2021", format = "%d-%m-%Y")) 


# we fix the rate of progression and rate of recovery 
# we estimate beta (transmission rate) and rho (reporting probability)

sigma = 1/5.1   # rate of progression
gamma = 1/2.1     # rate of recovery 

# population and seropositive populations are known
# we must estimate the initial number of infections 

n_pop = 15810388 # population of Guateng 
n_recov = round(n_pop * 0.562) # 56.2% seropositive  


# parameters to estimate are therfore:

pars = c("beta", "rho", "I0")

# modify data into a list format for Rstan ------------------------------------# 


stan_data_m1_solv = list(n_days = dim(sim_data)[1], # Number of observations 
                         n_pop = n_pop,
                         n_recov = n_recov,
                         y = sim_data$rep_inc[seed_omicron:dim(sim_data)[1]],
                         n_data = length(sim_data$rep_inc[seed_omicron:dim(sim_data)[1]]),
                         sigma = sigma,
                         gamma = gamma,
                         ts = sim_data$time, 
                         time_seed_omicron = seed_omicron)



# compile our Rstan model -----------------------------------------------------# 

m1_solv = stan_model(paste0(file_path,"models/model1_solverV1.stan"))


# run our Rstan model ---------------------------------------------------------# 


# randomly generate 3 different initial values from pluasible parameter range 
list_of_inits = list(
  ini_1(seed = 21),
  ini_1(seed = 05),
  ini_1(seed = 22)
)

# as model complexity increases, its helpful to keep track of the model run time  
time.start = Sys.time() 

m1_fit_solv = sampling(
  m1_solv,
  data = stan_data_m1_solv,
  chains=3,
  warmup=1000,
  iter=2000,
  init = ini_1,
  seed = 200522
)
time.end = Sys.time()

# model run time 

time.end - time.start 


# model diagnostics -----------------------------------------------------------# 

# We obtained no warnings about our model, great news! 

# To be on the safe side, lets check with some model diagnostics

check_divergences(m1_fit_solv) # no divergent transitions 

m1_solv_post= as.array(m1_fit_solv)


# plot markov chain trace plots to check for model convergence 

mcmc_trace(m1_solv_post, pars="lp__")
mcmc_trace(m1_solv_post, pars=pars)
mcmc_trace(m1_solv_post, pars="R_0")


# Univariate and bivariate marginal posterior distributions

pairs(m1_fit_solv, pars = pars, cex.labels=1.5, font.labels=9, condition = "accept_stat__")  


# overlay densisty estimates obtained for each chain separately 
mcmc_dens_overlay(m1_solv_post, pars=pars)


# summary  of the parameters we are interested in -----------------------------# 

m1_fit_solv_summary = summary(m1_fit_solv, pars = pars)$summary
print(m1_fit_solv_summary,scientific=FALSE,digits=2)




# Plotting the model output against the data ----------------------------------#     

m1_fit_solv_post = rstan::extract(m1_fit_solv)

plot_m1_solv_results  = m1_fit_solv_post$lambda %>% as.data.frame.table() %>%
  rename(ni = iterations, day = Var2, value = Freq) %>%
  dplyr::mutate(
    ni = as.numeric(ni),
    day = as.numeric(day)) %>% 
  group_by(day) %>% 
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  )  %>%  
  mutate(time = seed_omicron: dim(sim_data)[1]) %>% 
  left_join(sim_data) %>%  
  ggplot(aes(x = time , y = rep_inc)) +
  geom_point()+
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin=lower,ymax=upper))




