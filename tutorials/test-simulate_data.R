# install.packages("testthat") 

library(testthat)
context("Simulate data tetsing")
source("R/simulate_data.R")

test_that("simulate single var", {
  
 ODE_single_var = simulate_data_single_var(
    n_pop = 15810388 ,           
    immunity = 0.1,
    n_inf = 1 ,
    gamma = 1/4.17  ,              
    sigma = 1/3.03,   
    rho= 0.2,
    beta = 4,  
    ts= 1:100
  )
  
  # Test that the result is length 2 
  expect_that( dim(ODE_single_var), equals(c(100, 6)) );
  
  # Test that the first object returned is a data frame
  expect_that( is.data.frame(ODE_single_var[1]), equals(TRUE) );
 
  # Test that all values are positive 
  expect_that( all(ODE_single_var >= 0), equals(TRUE) )
  
})




test_that("simulate multi var", {
  
  ODE_multi_var = simulate_data_multi_var(
    n_pop = 15810388,               
    immunity = 0.5 , 
    n_inf_D =1,
    n_inf_O =1,
    gamma_D = 0.2,                
    sigma_D = 0.1,   
    gamma_O = 0.4,               
    sigma_O = 0.1, 
    rho_D  = 0.4,
    rho_O = 0.4,
    beta_D =1,
    beta_O = 2,
    epsilon = 0.1, 
    omega = 0.8,
    nu = 0.001, 
    ts = 1:100,
    time_int_start = 34,
    time_int_end = 88,
    time_seed_O = 12
  )
  
  # Test that the result is length 2 
  expect_that( dim(ODE_multi_var), equals(c(100, 10)) );
  
  # Test that the first object returned is a data frame
  expect_that( is.data.frame(ODE_multi_var[1]), equals(TRUE) );
  
  # Test that all values are positive 
  expect_that( all(ODE_multi_var >= 0), equals(TRUE) )
  
})


