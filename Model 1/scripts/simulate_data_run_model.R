library(deSolve)
library(dplyr)
library(rstan)
library(bayesplot)
library(loo)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

# generate data 
require(deSolve)

## First model is SEIR model 

SEIR <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+I+E+R
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma * E
    dI <- sigma * E - gamma*I 
    dR <- gamma*I
    
    
    return(list(c(dS,dE, dI, dR)))
  })
}

#### We assume an R0 of 9.5 - https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taac037/6545354?login=true
# We assume an incubation period of 5.1 days and an infectious period of 2.1 days

R0 = 9.5

# Sigma, the rate of progression from incubation to infectious is tehrefore 1/5.1 days 
sigma = 1 / 5.1

# gamma, the recovery rate is 1/2.1 days 
gamma = 1/2.1

# Beta, the transmission rate is therefore R0 * gamma 
beta = R0 * gamma 

params <- c(beta, sigma, gamma)


# Define our initial states
N = 15810388 # pop of Guateng 

# Seroprevalance survery in Guatang prior to 4th wave reported 56.2% immunity - https://www.medrxiv.org/content/10.1101/2021.12.20.21268096v1

R = round(N * 0.562)


## Assume 100 initial infections 
I = 100

initial_state <- c(S= N-R - I  , E = 0 , I=I , R=R)

times <- (1:153) # run in days for Sep-Jan 


# Solve the model using deSolve 
model <- ode(initial_state, times, SEIR, params)

out.df<-round(as.data.frame(model))

### Incidence is the rate that individuals enter the Infectious compartment 

# Inc(t) = sigma * E(t)
Inc = data.frame(Time = times,
                 Inc = round(out.df$E * sigma))

## plot epidemic curve 

ggplot(Inc, aes(x= Time , y = Inc)) +
  geom_point()


############  Fit model using Rstan  ############ 

# As the solve in Rstan can be very slow, we will instead estimate the changes in state at each time point at a fine temporal resolution

# As our model in days, we will update our estimate every 2.4 hours 

scale_time_step = 10


# Rstan requires the data in a list format 
data_rstan = list(n_obs = length(times), # Number of observations 
             n_pop = N,
             n_recov = R,
             y = Inc$Inc,
             scale_time_step = scale_time_step,
             sigma = sigma/scale_time_step,
             gamma = gamma/scale_time_step)



## Here we compile our Rstan model 

m1 <- stan_model( "Omicron_model_1.stan")

## These parameters tell Rstan how many chains we want to run, the number of iterations and how many to discard as burn in
n_chains=3
n_warmups=500
n_iter=5000

set.seed(1234)

# It helps the model to set initial values in a plausible range for each parameter 

# To start, we are going to fix the incubation and infectious period and estimate the initial number of cases and beta

# We think the R0 of Omicron could be between 5.5-12, so plausible starting values for Beta are 2-6

# 1-1000 seems like a reasonable range for the seed 

ini_1 = function(){
  list(beta=runif(1,2,12),
       I0 = runif(1,1,1000))}



time.start <- Sys.time()
model_fit_1 = sampling(
  m1,
  data = data_rstan,
  init = ini_1,
  chains = n_chains,
  warmup = n_warmups,
  iter = n_iter,
  seed = 13219
)
time.end <- Sys.time()
time.end - time.start 

model_fit_1_summary <- summary(model_fit_1, pars = c("lp__", "beta", "I0", "R_0"))$summary
print(model_fit_1_summary,scientific=FALSE,digits=2)

############## WAIC / LOO  ############## 

loo1 <- loo(model_fit_1, save_psis = TRUE)

