# Bayesian inference for SARS-CoV-2 transmission modelling, using Stan

### Author: Bethan Cracknell Daniels

The aim of this exemplar is to demonstrate how to design and fit a mathematical model of disease transmission to real data, in order to estimate key epidemiological parameters and inform public health responses. Specifically, we will model the emergence of the SARS-CoV-2 variant of concern Omicron in Gauteng, South Africa. To fit the model, we use Stan, a free, accessible and efficient Bayesian inference software. Adopting a Bayesian approach to model fitting allows us to account for uncertainty, which is especially important when modelling a new pathogen or variant. The transmission model uses compartments to track the populations movement between states, for instance from susceptible to infectious. By fitting a compartmental model to genomic and epidemiological surveillance data, we will recreate the transmission dynamics of Omicron and other circulating variants, and estimate key epidemiological parameters. Together these estimates are useful for guiding policy, especially in the early stages of an emerging variant or pathogen, when there are lots of unknowns.

## Prerequisites:

Required:
- Experience using R, for instance the Graduate school course *R programming*.
- Knowledge of Bayesian statistics, for instance, the textbook *A students guide to Bayesian statistics* by Ben Lambert is a great place to start. Chapter 16 also introduces Stan. 
- Download Rstan following [these](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) instructions. 

Beneficial:
- Some familiarity with Stan. 
- Experience of infectious disease modelling. 

## Useful external resources:

- [Bayesian workflow for disease transmission modeling in Stan](https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html#2_using_simulated_data_to_understand_our_model)
- [A students guide to Bayesian statistics is accompanied by a lecture course on youtube](https://www.youtube.com/playlist?list=PLwJRxp3blEvZ8AKMXOy0fc0cqT61GsKCG)
- [A brief introduction to Stan](https://www.alexpghayes.com/post/2018-12-24_some-things-ive-learned-about-stan/)
- [The Stan mannual](https://mc-stan.org/)
- [Statistical Rethinking](https://github.com/rmcelreath/stat_rethinking_2022) - this a thorough course on Bayesian data analysis which uses Stan. The course consists of lectures, homework and can be completed alongside reading the textbook of the same title

## Intended learning outcomes 

Upon completion of this tutorial, students will be able to:

1.	design an infectious disease compartmental model to answer public health questions. 
2.  compare methods of solving ODE using Stan. 
3.	write a Stan model to fit an infectious disease model.
4.  interpret Stan model diagnostics and implement appropriate solutions. 
5.  structure R code into files based on functionality. 
6.  write tests in R to check code. 

## Project structure 

This project is split into 3 chapters, contained in the folder *tutorials*:

- Chapter 1: Designing an infectious disease model
- Chapter 2: Fitting a single variant model in Stan
- Chapter 3: Fitting a multivariant model in Stan 


Each chapter is interactive, with questions and activities along the way. Therefore, for each chapter you will find 2 RMD files, one with and one without solutions. 

### Chapter 1 

This chapter will demonstrate how to design an infectious disease compartmental model, including choosing the model priors and likelihood. 


### Chapter 2 

This chapter will fit a single-variant compartmental model to simulated data in order to explore the transmission dynamics of omicron in Gauteng. In order to do fit the model, this chapter demonstrates how to: 

 1. simulate data 
 2. compare different methods of solving ordinary differential equations in Stan 
 3. code up an infectious disease model in Stan 
 4. fit infectious disease models in Stan 
 5. run model diagnostics 
 6. plot the model fit against data 
 
### Chapter 3
 
This chapter builds on everything introduced in chapter 2 in order to fit a more complicated model. Specifically, chapter 3 shows how to fit a multi-variant compartmental model to simulated data in order to explore the transmission dynamics of omicron and delta in Gauteng. The model accounts for population testing, vaccination and waning immunity. 
 
## Repository structure 

Each RMD file takes as in input multiple functions, each with a specific purpose. Frequently, these functions take as input the output of the function before. The functions are all contained in the folder *R*. There is also a *models* folder, which contains all compartmental models used in chapters 2 and 3. 

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


##### 

```diff 

@@ How to Use the Code @@

! Requirements

- R (version XXX or higher)
- RStudio (if needed)
- packages
  - rstan
  - others (for example ....)

! example:
+ XXXXXX
```
