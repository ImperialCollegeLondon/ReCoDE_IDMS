# Bayesian inference for SARS-CoV-2 transmission modelling, using Stan

<!-- buttons -->

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)

### Author: Bethan Cracknell Daniels

## Description

The aim of this exemplar is to demonstrate how to design and fit a mathematical model of disease transmission to real data, in order to estimate key epidemiological parameters and inform public health responses. Specifically, we will model the emergence of the SARS-CoV-2 variant of concern Omicron in Gauteng, South Africa. To fit the model, we use Stan, a free, accessible and efficient Bayesian inference software. Adopting a Bayesian approach to model fitting allows us to account for uncertainty, which is especially important when modelling a new pathogen or variant. The transmission model uses compartments to track the populations movement between states, for instance from susceptible to infectious. By fitting a compartmental model to epidemiological surveillance data, we will recreate the transmission dynamics of Omicron and other circulating variants, and estimate key epidemiological parameters. Together these estimates are useful for guiding policy, especially in the early stages of an emerging variant or pathogen, when there are lots of unknowns.

### Useful external resources:

- [Bayesian workflow for disease transmission modeling in Stan](https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html#2_using_simulated_data_to_understand_our_model)
- [A students guide to Bayesian statistics is accompanied by a lecture course on youtube](https://www.youtube.com/playlist?list=PLwJRxp3blEvZ8AKMXOy0fc0cqT61GsKCG)
- [A brief introduction to Stan](https://www.alexpghayes.com/post/2018-12-24_some-things-ive-learned-about-stan/)
- [The Stan manual](https://mc-stan.org/)
- [Statistical Rethinking](https://github.com/rmcelreath/stat_rethinking_2022) - this a thorough course on Bayesian data analysis which uses Stan. The course consists of lectures, homework and can be completed alongside reading the textbook of the same title

## Requirements:

### Academic

Required:
- Experience using R, for instance the Graduate school course *R programming*.
- Knowledge of Bayesian statistics, for instance, the textbook *A students guide to Bayesian statistics* by Ben Lambert is a great place to start. Chapter 16 also introduces Stan. 
- Download Rstan following [these](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) instructions. 

Beneficial:
- Some familiarity with Stan. 
- Experience of infectious disease modelling. 

### System

| Program                  | Version                  |
| ------------------------ | ------------------------ |
| R                        | Any                      |


## Learning outcomes 

Upon completion of this tutorial, students will be able to:

1.  design an infectious disease compartmental model to answer public health questions. 
2.  compare methods of solving ODE using Stan. 
3.  write a Stan model to fit an infectious disease model.
4.  interpret Stan model diagnostics and implement appropriate solutions. 
5.  structure R code into files based on functionality. 
6.  write tests in R to check code. 

## Getting Started
### Try the code with Docker container:

If you have Docker engine installed on your computer, please set Docker engine to use at least 6GB RAM, 4GB Swap, and 4 CPUs.

Pull the container image in a command window (Windows) or a terminal window (Mac or Linux) by running 

`docker pull jianlianggao/recode_idms:20220726`

When the image is pulled, to run the image in a container instance, please run the following command

`docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true jianlianggao/recode_idms:20220726`

Let the above command run in the terminal window and keep it open in the background, you can open `127.0.0.1:8787:8787` in your favourite web brower and you will have a RStudio (it is officially called `posit` now) interface to open and run tutorial code of chapter 1, chapter 2 or chapter 3 in R Markdown.  
You can modify the code in R Markdown but you do not have permissions to save the file. If you want to save changes you have made, please stop the container instance in the terminal window by pressing `Control + C`. Then start a new container instance by running

`docker run --rm -p 127.0.0.1:8787:8787 -v /tmp:/home/rstudio/data -e DISABLE_AUTH=true jianlianggao/recode_idms:20220703`

Again, please keep the command window (or terminal window) opened in the background.
View RStudio from your favourite web browser by visiting 127.0.0.1:8787 and now you should be able to save changes in the /home/rstudio/data folder, which is mapping to `/tmp` in your computer. You can copy the save files to other folder later on.

For Windows users, the mount path format `/tmp` is different, you may need to replace it with `d:/tmp` for example. This has been tested yet. We will update when we have a chance to test it on a Windows computer.
In RStudio in your web browser, please run config.R before running any other code.

## Project structure 

This project is split into 3 chapters:

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
 
This chapter builds on everything introduced in chapter 2 in order to fit a more complicated model. Specifically, chapter 3 shows how to fit a multistain compartmental model to simulated data in order to explore the transmission dynamics of omicron and delta in Gauteng. The model accounts for population testing, vaccination and waning immunity. Chapter 3 also includes some optional, open ended challenges. 

#### Repository Structure
The documents related to each chapter are organised in the following folders:

##### docs
This folder contains all the .md files associated with each chapter and was used to render the documentation.

The MD file testing introduces formal testing and the package *test that*. All *R.* files with the prefix *test* are scripts which run formal testing on a specific function. 

##### Tutorials

This folder contains the corresponding .html and .Rmd files of each chapter. 

##### R

Each RMD file takes as in input multiple functions, each with a specific purpose. Frequently, these functions take as input the output of the function before. This folder contains all the functions needed to run the Rmd file. All *R.* files with the prefix *test* are scripts which run formal testing on a specific function. 

As stated above, each function has a specific purpose and the functions are designed to be run in order. The functions are as follows: 

- *simulate_data.R*: Functions to produce simulated reported incidence data using the model *model1_deSolve.R* or *model2_deSolve.R*. Takes as input parameter values for the model. Outputs a data frame of solutions to the derivatives of all compartments at each time step. 

- *calc_sim_incidence.R*: Functions to calculate the simulated reported incidence. Takes as input a data frame of solutions to the derivatives of all compartments at each time step. Outputs a data frame and ggplot of reported incidence over time. 

- *run_stan_models.R*: Function to fit a stan model. Uses the function *draw_init_values.R*. At minimum, takes as in put a list of data to fit the model and a Stan model. Outputs a fitted stan model. 

- *draw_init_values.R*: Functions sourced within *run_stan_models.R*  which generates a different starting values for each Markov chain. Takes as input a seed value and the number of variants the model is fitting to. Outputs an initial value for each parameter and each chain. 

- *diagnose_stan_fit.R*: a function to run diagnostics on a Stan fit. Takes as input a fitted  Stan model and the parameters to check. Outputs the number of divergent transitions, diagnostic plots and parameter summary statistics. 

- *plot_model_fit.R*: a function to plot the results of a fitted Stan model against the data to which it was fit. Takes as input a fitted Stan model, the name of the variables to be plotted, and the simulated or observed data.

- *compare_param_est.R*: a function to compare  parameter estimates between models or methods of solving ODEs, i.e., in order to check whether a Stan model is able to recover true parameter estimates from simulated data. Takes as input a vector of true parameter values, estimated posterior mean and 95% CrI from a Stan fit and parameter names. Outputs plots comparing parameter values.  

##### models:

This folder contains the compartmental models used in part 1. 

**R models:**

- *model1_deSolve.R*: a function to solve a single strain SEIQR model which is sourced by the function *simulate_data.R*.

- *model2_deSolve.R*: a function to solve a multistrain SEIQRS model which is sourced by the function *simulate_data.R*. 

**Stan models (written in C++):**

- *model1_Euler_V1.stan*: a Stan model of a single strain SEIQR model, the ODEs are solved using the Euler method at a single day step. 

- *model1_Euler_V2.stan*: a Stan model of a single strain SEIQR model, the ODEs are solved using the Euler method at a user-defined time step.  

- *model1_RK_V1.stan*: a Stan model of a single strain SEIQR model, the ODEs are solved using the Runge-Kutta Method. 

- *model2_Euler_V1.stan*: a Stan model of a multistrain SEIQRS model, the ODEs are solved using the Euler method at a user-defined time step. 

- *model2_Euler_V2.stan*: a Stan model of a multistrain SEIQRS model, the ODEs are solved using the Euler method at a user-defined time step. 
