---
title: "Part 1"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_collapse: true

---





```r
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())
  
# Functions
source("R/simulate_data.R")
source("R/calc_sim_incidence.R")
source("R/draw_init_values.R")
source("R/run_stan_models.R")
source("R/diagnose_stan_fit.R")
source("R/plot_model_fit.R")
source("R/compare_param_est.R")
source("R/tidy_data.R")

# Models
m2_EU1 = stan_model("models/model2_Euler_V1.stan")
#m1_EU4 = stan_model("models/model1_Euler_V4.stan")
#m1_EU5 = stan_model("models/model1_Euler_V5.stan")
```


## Including Delta Variant

We are going to fit to both the third wave and fourth wave in Gauteng. We will assume that between April 2021 - October 2021 the incidence data is Delta, whilst between October and February the incidence data is Omicron. As before, we will seed both variants 1 month before fitting to the data. That means we seed Delta in March and we seed Omicron in September. 

During the third wave, driven by Delta, interventions were introduced between 15-06-21 and 13-09-2021. We will include these in our model. 

Defining key dates: 

```r
# Date we start fitting the incidence to Delta 
date_fit_D  =   as.Date.character("01-03-2021", format = "%d-%m-%Y")
 
# Date we start fitting the incidence data to Omicron 
date_fit_O = as.Date.character("01-10-2021", format = "%d-%m-%Y")

# Date we seed Delta (1 month before we fit)
date_seed_D =  as.Date.character("01-02-2021", format = "%d-%m-%Y") 

# Date we seed D (1 month before we fit)
date_seed_O = as.Date.character("01-09-2021", format = "%d-%m-%Y")

# Date interventions 
date_int= as.Date.character(c("15-06-2021","13-09-2021"), format = "%d-%m-%Y")

# Final data of modelling period 
end_date = as.Date.character("23-02-2022", format = "%d-%m-%Y")

# All dates of modelling period
all_dates_mv = seq.Date(from = date_seed_D, to = end_date ,  by = "days") # model times
```

Simulating multivariate data: 


We need redefine some of the variables and add some new variables. 

First, we redefine the number of recovered individuals, as we are now looking prior to the third wave. A sero-suvery in Gauteng reported a seroprevalence of 19% in at the end of January 2021 [5]. 

We will assume a single Omicron and Delta infection seeds the outbreaks and that the $R_0$ of both variants is 5. The probability of testing $\rho$ is the same as before, as are $\sigma$ and $\gamma$. 

We are now going to include vaccination in our model. The simplest way to model this is by assuming susceptible individuals enter the R compartment directly at a rate $\nu$.We will assume a constant per capita rate of vaccination of $\nu = 0.001$. 

We are also now going to account for waning immunity. Specifically, we will assume that, on average, immunity induced by an infection with Delta or through double vaccination lasts 90 days. We model this by including a new compartment called $SO$. Individuals enter this compartment from $R$ are a rate $$\epsilon = \frac{1}{180 days}$$  . From here, they be re-infected by Omicron, but not Delta. 

Let $\omega$ be the percentage reduction in transmission due to interventions introduced in the third wave during the dates defined above. We model this impact on transmission as follows: 


$$\beta' = \beta * (1-\omega)$$
We will assume that $\omega = 0.2$.  




Now, lets simulate our data: 

```r
# Define  variables 
n_pop = 15810388                  # population
  R0_D = 5
  R0_O = 6
  immunity_mv = 0.19 
  rho_D = 0.2
  rho_O = 0.05
  n_inf_D = 1
  n_inf_O = 1
  epsilon = 1/180
  omega = 0
  nu = 0.001
  gamma = 1/4.17                    # recovery rate 
  sigma = 1/3.03                    # latent rate 
  betaD = (R0_D * gamma) /(1-rho_D)
  betaO = (R0_O * gamma) /(1-rho_O)
```


```r
multi_var_sim_data = simulate_data_multi_var(
  immunity = immunity_mv, 
  n_inf_D = n_inf_D,
  n_inf_O = n_inf_O ,
  rho_D = rho_D,
  rho_O = rho_O,
  betaD =betaD,
  betaO = betaO,
  epsilon = epsilon, 
  omega = omega,
  nu = nu, 
  ts = 1:length(all_dates_mv),
  time_int_start = which(all_dates_mv == date_int[1]),
  time_int_end = which(all_dates_mv == date_int[2]),
  time_seed_O = which(all_dates_mv == date_seed_O)
)
```

```
## Warning: package 'deSolve' was built under R version 4.0.5
```

![](chapter-3_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-2.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-3.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-4.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-5.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-6.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-7.png)<!-- -->![](chapter-3_files/figure-html/unnamed-chunk-4-8.png)<!-- -->

And plot the reported incidence for Delta and Omicron: 


```r
multi_var_sim_inc = calc_sim_incidence_multi_var(
                             ODE_data = multi_var_sim_data,
                             all_dates = all_dates_mv,
                             date_fit_D = date_fit_D,
                             date_fit_O = date_fit_O,
                             rho_D = rho_D,
                             rho_O = rho_O)


multi_var_sim_inc[[2]] # Plot 
```

```
## Warning: Removed 28 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 242 row(s) containing missing values (geom_path).
```

![](chapter-3_files/figure-html/unnamed-chunk-5-1.png)<!-- -->





```r
# Removing missing data 
 y_D = multi_var_sim_inc[[1]]$rep_inc_D_noise[!is.na(multi_var_sim_inc[[1]]$rep_inc_D_noise)]
 y_O = multi_var_sim_inc[[1]]$rep_inc_O_noise[!is.na(multi_var_sim_inc[[1]]$rep_inc_O_noise)]

stan_fit_EU3_obs = run_stan_models(
  list_data =
    list(
      n_var = 2,                                                
      n_ts = length(all_dates_mv),                              # no. time steps 
      n_pop = n_pop,                                            # population  
      n_recov = round(immunity_mv*n_pop),                       # recovered population 
      I0 = c(n_inf_D, n_inf_O), 
      y_D = y_D,   # incidence data 
      y_O = y_O,   # incidence data
      n_data_D = length(y_D),  # no. data 
      n_data_O = length(y_O),  # no. data 
      sigma = sigma,                            # latent rate  
      gamma = gamma,                            # recovery rate 
      time_seed_O =                             # index seed omicron
        which(all_dates_mv == date_seed_O), 
      time_fit_D = 
         which(all_dates_mv == date_fit_D),
      time_fit_O = 
         which(all_dates_mv == date_fit_O),
      time_int_start = which(all_dates_mv == date_int[1]),
      time_int_end = which(all_dates_mv == date_int[2]),
      scale_time_step = 10,                                     # amount to reduce time step,
      nu  = nu
    ), 
  n_iter = 2000,
  n_warmup = 1000,
  model = m2_EU1,
  n_var = 2,
  model_no = 2
)
```

```
## Warning: The largest R-hat is NA, indicating chains have not mixed.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#r-hat
```

```
## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#bulk-ess
```

```
## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#tail-ess
```

```
## Time difference of 49.49258 mins
```


```r
m2_diag = diagnose_stan_fit(
  stan_fit_EU3_obs, #
  pars = c("beta[1]", "rho[1]","rho[2]", "beta[2]", "R_0[1]", "R_0[2]", "epsilon"))
```

```
## Warning: package 'bayesplot' was built under R version 4.0.5
```

```
## This is bayesplot version 1.8.1
```

```
## - Online documentation and vignettes at mc-stan.org/bayesplot
```

```
## - bayesplot theme set to bayesplot::theme_default()
```

```
##    * Does _not_ affect other ggplot2 plots
```

```
##    * See ?bayesplot_theme_set for details on theme setting
```

![](chapter-3_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
m2_diag
```

```
## $`markov chain trace plots`
```

![](chapter-3_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

```
## 
## $`univariate marginal posterior distributions`
```

![](chapter-3_files/figure-html/unnamed-chunk-7-3.png)<!-- -->

```
## 
## $`summary statistics of parameters`
##                mean      se_mean          sd        2.5%         25%
## beta[1]  1.22692344 0.0014021125 0.035974646 1.165455047 1.208020539
## rho[1]   0.20907376 0.0006676235 0.020378234 0.170725057 0.198934162
## rho[2]   0.03394768 0.0090271499 0.011599318 0.018485882 0.024464394
## beta[2]  2.62712953 0.8040353184 1.011017336 1.222158469 1.259454928
## R_0[1]   4.04360682 0.0041298083 0.025997236 3.993624028 4.029250539
## R_0[2]  10.62503585 3.3077950641 4.150889046 4.856257650 4.987686848
## epsilon  0.01328326 0.0046020502 0.006144651 0.005339107 0.005630429
##                 50%         75%       97.5%      n_eff     Rhat
## beta[1]  1.22006214  1.24212171  1.31509002 658.305295 1.021439
## rho[1]   0.20613294  0.21781727  0.25709108 931.686676 1.010579
## rho[2]   0.02942504  0.04824487  0.05177863   1.651063 3.105583
## beta[2]  3.09779248  3.41103762  3.89154287   1.581128 4.157505
## R_0[1]   4.03911886  4.05816751  4.10087105  39.627306 1.037122
## R_0[2]  12.60170974 13.84870640 15.74527220   1.574726 4.321027
## epsilon  0.01471992  0.01769424  0.02394656   1.782755 2.379327
```



```r
plot_model_fit_multi_var(
  stan_fit_EU3_obs,
  variable_model = "lambda_days",
  variable_data = "rep_inc_D_noise",
  data = multi_var_sim_inc[[1]])
```

```
## 
## Attaching package: 'tidyr'
```

```
## The following object is masked from 'package:rstan':
## 
##     extract
```

```
## `summarise()` has grouped output by 'time'. You can override using the `.groups`
## argument.Joining, by = c("time", "variant")
```

![](chapter-3_files/figure-html/unnamed-chunk-8-1.png)<!-- -->



```r
compare_param_est(
  parameter_names = c("beta[1]", "rho[1]","rho[2]", "beta[2]", "R_0[1]", "R_0[2]", "epsilon"),
  true_param_values = c(betaD, rho_D, rho_O, betaO , R0_D, R0_O, epsilon),
  param_values1 =  m2_diag[[3]][,c(1,4,8)],
  model_names = c("EU")
)
```

```
## [[1]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```
## 
## [[2]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

```
## 
## [[3]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-3.png)<!-- -->

```
## 
## [[4]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-4.png)<!-- -->

```
## 
## [[5]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-5.png)<!-- -->

```
## 
## [[6]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-6.png)<!-- -->

```
## 
## [[7]]
```

![](chapter-3_files/figure-html/unnamed-chunk-9-7.png)<!-- -->



# References 



- (1) Madhi SA, Kwatra G, Myers JE, et al. Population Immunity and Covid-19 Severity with Omicron Variant in South Africa. N Engl J Med 2022; 386(14): 1314-26.
- (2) Tanaka H, Ogata T, Shibata T, et al. Shorter Incubation Period among COVID-19 Cases with the BA.1 Omicron Variant. Int J Environ Res Public Health 2022; 19(10).
- (3) Lavezzo E, Franchin E, Ciavarella C, et al. Suppression of a SARS-CoV-2 outbreak in the Italian municipality of Vo’. Nature 2020; 584(7821): 425-9.
- (4) Liu, Y. & Rocklov, J. Liu Y, Rocklov J. The reproductive number of the Delta variant of SARS-CoV-2 is far higher compared to the ancestral SARS-CoV-2 virus. J Travel Med 2021; 28(7).
- (5) Mutevedzi PC, Kawonga M, Kwatra G, et al. Estimated SARS-CoV-2 infection rate and fatality risk in Gauteng Province, South Africa: a population-based seroepidemiological survey. Int J Epidemiol 2022; 51(2): 404-17.
