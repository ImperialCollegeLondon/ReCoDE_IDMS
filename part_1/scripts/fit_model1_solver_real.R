################################################################################
# This code fits model 1 to the observed reported incidence data for the       #
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

source(paste0(file_path, "R/model1_functions.R"))

# read in observed data to fit the model to -----------------------------------#

real_data = read.csv(paste0(file_path, "data/reported_cases.csv"))


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

sigma = 1 / 5.1   # rate of progression
gamma = 1 / 2.1     # rate of recovery

# population and seropositive populations are known
# we must estimate the initial number of infections

n_pop = 15810388 # population of Guateng
n_recov = round(n_pop * 0.562) # 56.2% seropositive


# parameters to estimate are therfore:

pars = c("beta", "rho", "I0")

# modify data into a list format for Rstan ------------------------------------#

real_data$date = as.Date.character(real_data$date,  format = "%Y-%m-%d")
y = real_data %>%
  filter(
    date >= as.Date.character("01-10-2021", format = "%d-%m-%Y") &
      date <=  as.Date.character("23-02-2022", format = "%d-%m-%Y")
  ) %>%
  select(cases) %>%
  mutate(time = c(seed_omicron:length(all_dates)))

plot(y$cases)

stan_inc_m1_solv = list(
  n_days = length(all_dates),
  n_pop = n_pop,
  n_recov = n_recov,
  y = y$cases,
  n_data = length(y$cases),
  sigma = sigma,
  gamma = gamma,
  ts = 1:length(all_dates),
  time_seed_omicron = seed_omicron
)



# compile our Rstan model V1 --------------------------------------------------#

# This is the same model we fit to our simulated date 

m1_solv = stan_model(paste0(file_path, "models/modelV1_solver.stan"))


# run our Rstan model V1 ------------------------------------------------------#


time.start = Sys.time()
m1_fit_solv_real = sampling(
  m1_solv,
  data = stan_inc_m1_solv,
  chains = 3,
  warmup = 1000,
  iter = 2000,
  init = ini_1,
  seed = 2105200
)
time.end = Sys.time()

# model run time

time.end - time.start

# Running our model with the real data produced a lot of warnings:

#-------------------------------------------------------------------------------
# Warning message:
#   File monitoring failed for project at "Q:/ReCoDE_IDMS"
# Error 2 (The system cannot find the file specified)
# Features disabled: R source file indexing, Diagnostics
# Warning messages:
#   1: There were 1373 transitions after warmup that exceeded the maximum treedepth.
# Increase max_treedepth above 10. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
# 2: Examine the pairs() plot to diagnose sampling problems
#
# 3: The largest R-hat is NA, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#r-hat
# 4: Bulk Effective Samples Size (ESS) is too low,
#indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess
# 5: Tail Effective Samples Size (ESS) is too low,
#indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess

#-------------------------------------------------------------------------------

# If you are not familiar with these warnings, I recommend following the links
# above to read more about what they mean and why they are worrying.


# Note too, our model run time increased from a few minutes to ~50 mins.

# Let's take a look at the model diagnostics again to work out what is going
# wrong!



# model diagnostics -----------------------------------------------------------#


# still no divergent transitions, that is good news.

check_divergences(m1_fit_solv_real)

m1_solv_post_real = as.array(m1_fit_solv_real)

# The largest R-hat was NA, indicating chains have not mixed, lets take a look.


mcmc_trace(m1_solv_post_real, pars = "lp__")
mcmc_trace(m1_solv_post_real, pars = pars)
mcmc_trace(m1_solv_post_real, pars = "R_0")


# Univariate and bivariate marginal posterior distributions

pairs(
  m1_fit_solv_real,
  pars = pars,
  cex.labels = 1.5,
  font.labels = 9,
  condition = "accept_stat__"
)


# overlay density estimates obtained for each chain separately
mcmc_dens_overlay(m1_fit_solv_real, pars = pars)

# -----------------------------------------------------------------------------#

# All these diagnostics indicate that we have a multimodal distribution, i.e.,
# there are at least two posterior modes.

# In particular, looking at the rho parameter (the probability of detecting and
# reporting a case), we can see that the model is exploring areas of both very
# high and very rho values.

# Currently we assume a beta(1,1) prior for rho, which is uniform between 0 and
# 1 (the required bounds of a probability parameter):

# we can check this by simulating 1000000 random samples from a beta
# distribution with shape1 = 1, shape2 = 1.

sample_b = rbeta(1000000, 1, 1)

plot(density(sample_b))

# However, a seroprevalence study in Guateng found that cases were under-reported
# by 2-20 fold: https://academic.oup.com/ije/article/51/2/404/6414575

# From this, we expect rho to most likely be between 0.05-0.5.

# This is prior knowledge not currently reflected in our model!


################################################################################
# Activity --------------------------------------------------------------------#
################################################################################

# (1) Save the existing stan model "model1_solverV1.stan" as 
#     "model1_solverV2.stan"

# (2) Change the beta prior on rho so that it reflects our prior knowledge that
#     we think rho is mostly likely less than 0.5. Hint: use the rbeta function
#     or the website https://ben18785.shinyapps.io/distribution-zoo/ to explore
#     the probability density function of the beta distribution using different
#     values of shape1 and shape2

# (3) Rerun the model using the code below.



# compile our second Rstan model ----------------------------------------------#

m1_solV2 = stan_model(paste0(file_path, "models/model1_solverV2.stan"))


# run our Rstan model ---------------------------------------------------------#


time.start = Sys.time()
m1_fit_solv_real2 = sampling(
  m1_solV2,
  data = stan_inc_m1_solv,
  chains = 3,
  warmup = 1000,
  iter = 2000,
  init = ini_1,
  seed = 2105200
)
time.end = Sys.time()

# model run time

time.end - time.start


# Plotting the model output against the data ----------------------------------#

m1_fit_solv_post_real = rstan::extract(m1_fit_solv_real)

plot_m1_solv_results_real  = m1_fit_solv_post_real$lambda %>% as.data.frame.table() %>%
  rename(ni = iterations, day = Var2, value = Freq) %>%
  dplyr::mutate(ni = as.numeric(ni),
                day = as.numeric(day)) %>%
  group_by(day) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  )  %>%
  mutate(time = seed_omicron:length(all_dates)) %>%
  left_join(y) %>%
  ggplot(aes(x = time , y = cases)) +
  geom_point() +
  geom_line(aes(y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper))




# summary  of the parameters we are interested in -----------------------------#

m1_fit_solv_summary = summary(m1_fit_solv, pars = pars)$summary
print(m1_fit_solv_summary,
      scientific = FALSE,
      digits = 2)
