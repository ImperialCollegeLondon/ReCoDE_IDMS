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

# model times 
ts = 1:length(all_dates)

# modify data into a list format for Rstan ------------------------------------# 


stan_data_m1_solv = list(n_days = length(times), # Number of observations 
                         n_pop = N,
                         n_recov = R,
                         y = Inc$Inc[seed_omicron:length(Inc$Inc)],
                         n_data = length(Inc$Inc[seed_omicron:length(Inc$Inc)]),
                         ts = times,
                         sigma = sigma,
                         gamma = gamma,
                         time_seed_omicron = seed_omicron)



## Here we compile our Rstan model 

m1_solv =stan_model( "models/Omicron_model_1_solver.stan")


time.start <- Sys.time()
m1_fit_solv = sampling(
  m1_solv,
  data = stan_data_m1_solv,
  init = ini_1,
  n_chains=3,
  n_warmups=1000,
  n_iter=2000,
  n_thin = 5
)
time.end <- Sys.time()
time.end - time.start 


## Summary of the parameters we are interested in 

m1_fit_solv_summary <- summary(m1_fit_solv, pars = c("lp__", "beta", "I0", "R_0", "rho"))$summary
print(m1_fit_solv_summary,scientific=FALSE,digits=2)

############## WAIC / LOO  ############## 

loom1_solv = loo(m1_fit_solv, save_psis = TRUE)



################################      Plotting the fits ################################     

m1_fit_solv_post = rstan::extract(m1_fit_solv)

m1_solv_results  = m1_fit_solv_post$lambda %>% as.data.frame.table() %>%
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
  mutate(Time = seed_omicron : length(times)) %>% 
  left_join(Inc) %>%  
  ggplot(aes(x = Times, y = Inc)) +
  geom_point()+
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin=lower,ymax=upper))



