library(deSolve)
library(dplyr)
library(rstan)
library(bayesplot)
library(loo)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

# generate data 
require(deSolve)


############  Dates of model ############ 

time_intervention = c("01-07-2021", "13-09-2021")
time_omicron_emerge = "01-08-2021"
time_seed_omicron = "01-09-2021"
time_seed_delta = "01-04-2021"
start_date = "01-03-2021"
end_date = "23-02-2022"

start_date_d = as.Date.character(start_date, format = "%d-%m-%Y")
end_date_d = as.Date.character(end_date, format = "%d-%m-%Y")

daily_date = seq.Date(from = start_date_d, to = end_date_d ,  by = "days")

time_switch1 =  which(daily_date == as.Date.character(time_intervention[1], format = "%d-%m-%Y")) 
time_switch2 = which(daily_date == as.Date.character(time_intervention[2], format = "%d-%m-%Y")) 
seed_omicron =  which(daily_date == as.Date.character(time_seed_omicron, format = "%d-%m-%Y"))  
emerge_omicron =  which(daily_date == as.Date.character(time_omicron_emerge, format = "%d-%m-%Y")) 
seed_delta =  which(daily_date == as.Date.character(time_seed_delta, format = "%d-%m-%Y"))  
times = 1:length(daily_date)

R0_D = 3
R0_O = 5

# We assume an incubation period of 5.1 days and an infectious period of 2.1 days for both variants 


# Sigma, the rate of progression from incubation to infectious is tehrefore 1/5.1 days 
sigma = 1 / 5.1

# gamma, the recovery rate is 1/2.1 days 
gamma = 1/2.1

# Assume a reporting rate of 30% and that interventions reduce transmission by 50%
rho = 0.3
omega = 0.5


# Beta, the transmission rate is therefore R0 * gamma 
betaD = (R0_D * gamma) / (1-rho)
betaO = (R0_O * gamma) / (1-rho)


# Define our initial states
n_pop = 15810388 # pop of Guateng 

# Seroprevalance survery in Guatang prior to 4th wave reported 56.2% immunity - https://www.medrxiv.org/content/10.1101/2021.12.20.21268096v1

n_recov = round(n_pop * 0.562)


## Assume 100 initial infections of both 
I0 = c(100,100)

initial_state2 = c(S= n_pop-n_recov -I0[1], ED = 0 , ID=I0[1] , EO = 0, IO = 0, Q=0, R=n_recov)



params2 = c(
  betaD = betaD,
  betaO = betaO,
  sigma = sigma,
  gamma = gamma,
  rho = rho,
  omega = omega,
  time_switch = time_switch1,
  time_switch_off =time_switch2
)




seed_omicron_event =  data.frame(var = c("IO", "S"),
                                 time = emerge_omicron,
                                 value = c(I0[2],-I0[2]),
                                 method = c("add" ,"add"))


## First model is SEIR model 


SEIQR <- function(times, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S + ID + ED + EO + IO + R + Q
    
    ##  model interventions
    
    
    if(times >= time_switch & times <= time_switch_off) foiD = (1-omega) * betaD*(ID/N)
    else foiD = betaD*(ID/N)
    
    
    if(times %in% time_switch:time_switch_off) foiO = (1-omega) * betaO*(IO/N)
    else foiO = betaO*(IO/N)
    
    
    dS = - (foiD+foiO)*S 
    dED = foiD * S- sigma * ED
    dEO = foiO * S - sigma * EO
    
    dID = (1-rho) * sigma * ED - gamma * ID 
    dIO = (1-rho) * sigma * EO - gamma * IO 
    
    dQ = rho * sigma * (EO + ED) - gamma * Q
    
    dR = gamma* (ID+ IO + Q)
    
    
    return(list(c(dS,dED, dID,dEO, dIO,dQ, dR )))
  })
}


############   Solve the model using deSolve ############  
model2 = data.frame(ode(initial_state2, times, SEIQR, params2, events = list(data = seed_omicron_event)))
print(model2,scientific=FALSE,digits=4)

ggplot(model2, aes(x= time , y = IO)) +
  geom_point()



### Incidence is the rate that individuals enter the Infectious compartment 

# Inc(t) = sigma * E(t)
Inc = data.frame(Time = times,
                 IncD = round((1-rho) * model2$ED * sigma),
                 IncO = round((1-rho) * model2$EO * sigma))

## plot epidemic curve 

ggplot(Inc, aes(x= Time , y = IncO)) +
  geom_point()


############  Fit model using Rstan  ############ 

# As the solve in Rstan can be very slow, we will instead estimate the changes in state at each time point at a fine temporal resolution

# As our model in days, we will update our estimate every 2.4 hours 

scale_time_step = 10


# Rstan requires the data in a list format 
data_rstan2 = list(n_days = length(times), # Number of observations 
                   n_var = 2, # Number of variants 
                   n_pop = n_pop,
                   n_recov = n_recov,
                   n_ts = length(times) * scale_time_step, 
                   n_dataD = length(Inc$IncD[-c(1:seed_delta)]),
                   n_dataO = length(Inc$IncO[-c(1:seed_omicron)]), 
                   y_D = Inc$IncD[-c(1:seed_delta)]   ,
                   y_O = Inc$IncO[-c(1:seed_omicron)],
                   scale_time_step = scale_time_step,
                   sigma = sigma/scale_time_step,
                   gamma = gamma/scale_time_step,
                   gamma_days = gamma,
                   omega = omega,
                   time_switch = time_switch1 * scale_time_step,
                   time_switch_off =time_switch2 * scale_time_step,
                   time_seed_omicron = seed_omicron * scale_time_step)



## Here we compile our Rstan model 

m2 = stan_model( "Omicron_model_2.stan")

## These parameters tell Rstan how many chains we want to run, the number of iterations and how many to discard as burn in
n_chains=3
n_warmups=500
n_iter=1000

set.seed(1234)

# It helps the model to set initial values in a plausible range for each parameter 

# To start, we are going to fix the incubation and infectious period and estimate the initial number of cases and beta

# We know that the outbreaks of delta and omicron take off, so an R0 between 1 and 8 seems reasoable
# This translates to starting values of beta betwen 0.5-4

# 1-1000 seems like a reasonable range for the seed 

ini_2 = function(){
  list(beta=runif(2,1.5,3),
       I0 = runif(2,1,100),
       rho = runif(1,0.2,0.4), 
       omega = runif(1,0.4,0.6))}



time.start <- Sys.time()
model_fit_2= sampling(
  m2,
  data = data_rstan2,
  init = ini_2,
  chains = n_chains,
  warmup = n_warmups,
  iter = n_iter,
  seed = 13219,
  control = list(
    adapt_delta = 0.9, 
    max_treedepth = 15
  )
)
time.end <- Sys.time()
time.end - time.start 

model_fit_2_summary <- summary(model_fit_2, pars = c("lp__", "rho", "omega", "I0[2]", "I0[1]", "R_0[1]", "R_0[2]"))$summary
print(model_fit_2_summary,scientific=FALSE,digits=2)


print(summary(model_fit_2, pars = c("lp__", "I"))$summary)

############## WAIC / LOO  ############## 

loo2 <- loo(model_fit_2, save_psis = TRUE)


############## plot model output vs. data ############## 

model2_fit <- rstan::extract(model_fit_2)

sim_data = data.frame(Y = c(Inc$IncD,Inc$IncO),
                      variant = rep(c("Delta", "Omicron"), each = n_obs))

pred_df = model2_fit$lambda_days %>% as.data.frame.table() %>%  
  rename(np = iterations, day = Var2, variant = Var3, value = Freq) %>%  
  dplyr::mutate(
    np = as.numeric(np),
    day = as.numeric(day),
    variant = factor(variant, labels = c("Delta", "Omicron"))) %>% 
  group_by(day, variant) %>% 
  summarise(
    lower = quantile(value, 0.025, na.rm =T),
    mean = mean(value, na.rm =T),
    upper = quantile(value, 0.975, na.rm =T)
  ) 


df = cbind(data_rstan2$y_D, data_rstan2$y_O) %>% as_data_frame() %>%
  rename("Delta" = y_D, "Omicron" = y_O) %>% 
  tidyr::gather(variant, value) %>%
  mutate(variant = factor(variant, levels = c("Delta", "Omicron"))) %>% 
  group_by(variant) %>% 
  mutate(day = row_number()) %>%
  ungroup() 


ggplot(data = pred_df) +
  geom_ribbon(aes(x = day, ymin = lower, ymax = upper, fill = variant), alpha = 0.3) +
  geom_line(aes(x = day, y = mean, colour = variant)) +
  geom_point(data = df, aes(x = day, y = value, colour = variant)) +
  ylab("Reported Incidence") + xlab("Day")
