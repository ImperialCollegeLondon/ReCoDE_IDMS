# Bayesian inference for SARS-CoV-2 transmission modelling, using Stan: Part 1 

## Author: Bethan Cracknell Daniels


In part 1, you will find 2 RMD files:

- ***part1.Rmd:*** this is the main tutorial which talks through: 

 1. designing an infectious disease model
 2. coding it up in Stan 
 3. simulating data 
 4. fitting the model to data
 5. model diagnostics
 6. plotting the model fit
 
 
Broadly, each one of these activities is achieved by one or more functions, which mostly takes as input the output of the function before. The Rmd. file is also designed to be interactive, with questions and activities along the way.

- ***part1_Solutions.Rmd:*** is as above, but with the solutions. 



Also in the part 1 you will find the following folders: 

### R:

This folder contains all the functions needed to run the Rmd file. As stated above, each function has a specific purpose and the functions are designed to be run in order. The functions are as follows: 

- *simulate_data.R*: A function to produce simulated reported incidence data using the model *model1_deSolve.R*. Takes as input parameter values for the model. Outputs a data frame of solutions to the derivatives of all compartments at each time step. 

- *calc_sim_incidence.R*: A function to calculate the simulated reported incidence. Takes as input a data frame of solutions to the derivatives of all compartments at each time step. Outputs a data frame and ggplot of reported incidence over time. 

- *run_stan_models.R*: A function to fit a stan model. Uses the function *draw_init_values.R*. At minimum, takes as in put a list of data to fit the model and a stan model. Outputs a fitted stan model. 

- *draw_init_values.R*: A function sourced within *run_stan_models.R*  which generates a different starting values for each Markov chain. Takes as input a seed value and outputs an initial value for each parameter and each chain. 

- *diagnose_stan_fit.R*: a function to run diagnostics on a Stan fit. Takes as input a fitted  stan model and the parameters to check. Outputs the number of divergent transitions, diagnostic plots and parameter summary statistics. 

- *plot_model_fit.R*: a function to to plot the results of a fitted stan model against the data to which it was fit. Takes as input a fitted Stan model, the name of the variables to be plotted, and the simulated or observed data.

### models:

This folder contains the compartmental models used in part 1. 

**R models:**

- *model1_deSolve.R*: a function to solve SEIQR model which is sourced by the function *simulate_data.R*. 

**Stan models (written in C++):**

- *model1_Euler_V1.stan*: a Stan model which solves the ODEs using the Euler method at a single day step. 

- *model1_Euler_V2.stan*: a Stan model which solves the ODEs using the Euler method at a used defined time step. 

- *model1_Euler_V3.stan* and *model1_Euler_V4.stan* are modified versions of *model1_Euler_V2.stan*.

- *model1_RK_V1.stan*: a Stan model which solves the ODEs using the Runge-Kutta Method. 


### images:


This file contains figures used in the RMD as learning resources, and can largely be ignored. 

