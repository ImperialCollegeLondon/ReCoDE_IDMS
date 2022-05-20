################################################################################
# This code simulates reported incidence data for the emergence of Omicron in  #
# Gauteng, South Africa. The code uses DeSolve to solve an SEIQR model,        #
# described in the README.md filer, under part 1.                              #
################################################################################



# set up ----------------------------------------------------------------------#  

rm(list = ls())
file_path = "part_1/"

# load packages and source user defined functions -----------------------------#       

library(deSolve)
library(tidyverse)

source(paste0(file_path,"R/model1_functions.R"))

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

# model parameters 

R0 = 5     # reproduction number 
sigma = 1/5.1   # rate of progression
gamma = 1/2.1     # rate of recovery 
rho = 0.2         # reporting probability 

beta =( R0 * gamma )/ (1-rho)    # beta, the transmission rate = R0 / (1-rho) * gamma 



# initial states 

n_pop = 15810388 # population of Guateng 
n_recov = round(n_pop * 0.562) # 56.2% seroprev 
n_inf = 10 # 1 initial infection 
S0 = n_pop - n_recov # susceptible at time of model fittin

initial_state = c(S= S0 - n_inf  , E = 0 , I=n_inf ,Q=0, R=n_recov)

params = c(beta, sigma, gamma,rho,S0)

# Solve the model using deSolve -----------------------------------------------#   

model = ode(initial_state, ts, SEIQR, params)

out.df  = round(as.data.frame(model))

# plot the model states to check they make sense ------------------------------#  

for(i in 2:ncol(out.df)) plot(out.df[,i])

# as an extra sanity check, lets run the model again with R0 = 1---------------#  

beta_check  =( 1 * gamma )/ (1-rho) 

params_check = c(beta=beta_check, sigma, gamma,rho)

model_check = ode(initial_state, ts, SEIQR, params_check)

# what do we expect to happen?

for(i in 2:ncol(model_check)) plot(model_check[,i])

# as expected, R0 = 1 means that, on average, each infectious person infects 
# just one more person. This means that the epidemic doesn't take off. 


# going back to our simulated epidemic, lets calculate the reported incidence   

# Inc(t) =  rho * sigma * E(t)
# (i.e., the rate of entry into the Q compartment)

sim_data = data.frame(time = ts,
                      rep_inc = round(rho *out.df$E * sigma))

# plot epidemic curve ---------------------------------------------------------#  

ggplot(sim_data, aes(x= time , y = rep_inc)) +
  geom_point()


# save simulated data ---------------------------------------------------------# 

write.csv(sim_data, paste0(file_path, "data/sim_data.csv"))

