# Bayesian inference for SARS-CoV-2 transmission modelling, using Stan

## Author: Bethan Cracknell Daniels

The aim of this exemplar is to demonstrate how to design and fit a mathematical model of disease transmission to real data, in order to estimate key epidemiological parameters and inform public health responses. Specifically, we will model the emergence of the SARS-CoV-2 variant of concern Omicron in Gauteng, South Africa. To fit the model, we use Stan, a free, accessible and efficient Bayesian inference software. Adopting a Bayesian approach to model fitting allows us to account for uncertainty, which is especially important when modelling a new pathogen or variant. The transmission model uses compartments to track the population’s movement between states, for instance from susceptible to infectious. By fitting a compartmental model to genomic and epidemiological surveillance data, we will recreate the transmission dynamics of Omicron and other circulating variants, and estimate key epidemiological parameters. Together these estimates are useful for guiding policy, especially in the early stages of an emerging variant or pathogen, when there are lots of unknowns

### Prerequisites:

Required:
•	Experience using R, for instance the Graduate school course *R programming *. 
•	Knowledge of Bayesian statistics, fr isntance, the textbook *A students guide to Bayesian statistics* by Ben Lambert is a great place to start. Chapter 16 also introduces Stan. 

Beneficial:
•	Some familiarity with Stan. 

### Project learning outcomes:

1.	Students will be able to design an infectious disease compartmental model to answer public health questions. 
2.	Students will gain appreciation of the benefits of using simulated data to debug models.  
3.	Students will be able to code up an Rstan model to fit an infectious disease model to observed data and estimate parameters.
4.	Students will be able to clean and format data using Tidyverse and produce high level plots using ggplot, in R. 

### Project structure 

The project is split into XX sections, with each section introducing a progressibely more comlex infectious disease model. Within each section there are the following modules: 

• Designing an infectious disease model 
•	Cleaning data in R
•	Simulating data in R
•	Model fitting in Rstan
•	Model diagnostics 
•	Plotting results in R

Excluding the first, each module will contain the following folders:

- data
- scripts
- R (functions)
- results / figures 

Moreover, each module will be accompanied  by lots of internal and external resources! 

