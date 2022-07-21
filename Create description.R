-library(usethis)
use_description(check_name = F,
                roxygen = F,
                fields = list(
                  Title = "Bayesian inference for SARS-CoV-2 transmission modelling, using Stan",
                  Description = "A project to fit a mathematical model of disease transmission to epidemiological data, using R and Stan. 
                  To learn more about this project, start with README.md",
                  `Authors@R` =  person("Bethan", "Cracknell Daniels", email = "bnc19@ic.ac.uk", 
                         role = c("aut", "cre")),
                  Version = "1.0"
               )
)
usethis::use_package("bayesplot", min_version ="1.8.1")
usethis::use_package("deSolve" , min_version = "1.2.8")
usethis::use_package("dplyr", min_version = "1.0.6")
usethis::use_package("ggplot2" , min_version = "3.3.6")
usethis::use_package("rstan" , min_version = "2.21.5")
usethis::use_package("tidyr" , min_version = "1.1.3")

usethis::use_package("R" , min_version = "4.0.4", "Depends")
