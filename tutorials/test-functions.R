# set up -----------------------------------------------------------------------

# install.packages("testthat") 
# setwd("C:/Users/bnc19/Desktop/ReCoDE_IDMS/tutorials")

# packages 
library(testthat)
context("Testing")

# source functions to test 
source("R/simulate_data.R")
source("R/draw_init_values.R")
source("R/compare_param_est.R")

# test functions from simulate_data.R ------------------------------------------
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
  
  # Test that the dimensions are correct 
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


# test draw initial values -----------------------------------------------------
test_that("draw init values 1", {
 
  n_var = 2 
  init_values = draw_init_values(seed = 10, n_var = n_var)

# Test that rho is between 0 and 1 
expect_that( all(init_values$rho > 0 & init_values$rho < 1), equals(TRUE) );

# Test that beta values > 0  
expect_that( all(init_values$beta > 0), equals(TRUE) );

# Test that kappa values > 0  
expect_that( all(init_values$k > 0), equals(TRUE) );

# Test that beta values returned equal number of variants 
expect_that( length(init_values$beta), equals(n_var) );

# Test that the first object returned is a list
expect_that( is.list(init_values), equals(TRUE) );


})

# ------------------------------------------------------------------------------


  test_that("compare_param_est", {
    
    param_values1 = data.frame(mean = c(1,2),
                               lower = c(0.4,0.6),
                               upper = c(0.8,0.9))
    
    comp = compare_param_est(
    parameter_names = c("rho", "beta"),
    true_param_values = c(1,3),
    param_values1 = param_values1,
    model_names = c("EU")
    )
    
    # Test that the object returned is a list
    expect_that( is.list(comp), equals(TRUE) );
    
    # Test that each element of list is a ggplot 
    for(i in 1:dim(param_values1)[1]){
    expect_that( is.ggplot(comp[[i]]), equals(TRUE) );
    }
    
  })
