# Bayesian inference for SARS-CoV-2 transmission modelling, using Stan

## Author: Bethan Cracknell Daniels

The aim of this exemplar is to demonstrate how to design and fit a mathematical model of disease transmission to real data, in order to estimate key epidemiological parameters and inform public health responses. Specifically, we will model the emergence of the SARS-CoV-2 variant of concern Omicron in Gauteng, South Africa. To fit the model, we use Stan, a free, accessible and efficient Bayesian inference software. Adopting a Bayesian approach to model fitting allows us to account for uncertainty, which is especially important when modelling a new pathogen or variant. The transmission model uses compartments to track the population’s movement between states, for instance from susceptible to infectious. By fitting a compartmental model to genomic and epidemiological surveillance data, we will recreate the transmission dynamics of Omicron and other circulating variants, and estimate key epidemiological parameters. Together these estimates are useful for guiding policy, especially in the early stages of an emerging variant or pathogen, when there are lots of unknowns

### Prerequisites:

Required:
- Experience using R, for instance the Graduate school course *R programming*.
- Knowledge of Bayesian statistics, for instance, the textbook *A students guide to Bayesian statistics* by Ben Lambert is a great place to start. Chapter 16 also introduces Stan. 
- Download Rstan following [these](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) instructions. 

Beneficial:
- Some familiarity with Stan. 
- Experience of infectious disease modelling. 

### Useful external resources:

- [Bayesian workflow for disease transmission modeling in Stan](https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html#2_using_simulated_data_to_understand_our_model)
- [A students guide to Bayesian statistics is accompanied by a lecture course on youtube](https://www.youtube.com/playlist?list=PLwJRxp3blEvZ8AKMXOy0fc0cqT61GsKCG)
- [A brief introduction to Stan](https://www.alexpghayes.com/post/2018-12-24_some-things-ive-learned-about-stan/)
- [The Stan mannual](https://mc-stan.org/)
- [Statistical Rethinking](https://github.com/rmcelreath/stat_rethinking_2022) - this a thorough course on Bayesian data analysis which uses Stan. The course consists of lectures, homework and can be completed alongside reading the textbook of the same title

### Intended learning outcomes 

Upon completion of this tutorial, students will be able to:

1.	design an infectious disease compartmental model to answer public health questions. 
2.  compare methods of solving ODE using Stan. 
3.	write a Stan model to fit an infectious disease model.
4.  interpret Stan model diagnostics and implement appropriate solutions. 
5.  structure R code into files based on functionality. 
6.  write tests in R to check code. 

### Project structure 

The project is split into 2 parts, with part 2 building on the model and fitting process introduced in part 1. Within each part, the tutorial will design an infectious disease model, code it up in Stan, fit the model to simulated and real data and run model diagnostics and plot model outputs. 


