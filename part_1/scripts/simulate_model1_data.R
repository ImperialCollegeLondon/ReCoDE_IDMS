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

R0 = 9.5          # reproduction number 
sigma = 1 / 5.1   # rate of progression
gamma = 1/2.1     # rate of recovery 
beta = R0 * gamma # beta, the transmission rate = R0 * gamma 
rho = 0.2         # reporting probability 

params = c(beta, sigma, gamma,rho)


# initial states 

n_pop = 15810388 # population of Guateng 
n_recov = round(n_pop * 0.562) # 56.2% seroprev 
n_inf = 1 # 1 initial infection 

initial_state = c(S= n_pop - n_recov - n_inf  , E = 0 , I=n_inf ,Q=0, R=n_recov)


# Solve the model using deSolve -----------------------------------------------#   

model = ode(initial_state, ts, SEIQR, params)

out.df  = round(as.data.frame(model))

# calculate the reported incidence --------------------------------------------#  

# Inc(t) =  rho * sigma * E(t)
# (i.e., the rate of entry into the Q compartment)

sim_data = data.frame(time = ts,
                      rep_inc = round(rho *out.df$E * sigma))

# plot epidemic curve ---------------------------------------------------------#  

ggplot(sim_data, aes(x= time , y = rep_inc)) +
  geom_point()


# save simulated data ---------------------------------------------------------# 

write.csv(sim_data, paste0(file_path, "data/sim_data.csv"))
