 
 // We first need to define a function which solves our compartmental model
 // and outputs the derivatives of all compartments at each time step 
 
 functions {
real[] SIR(real t,    // time
  real[] y,           // system state 
  real[] theta,       // parameters 
  real[] x_r,         // real data
  int [] x_i) {       // integer data
  
  
  // define integer data // 

  int n_recov = x_i[1]; 
  int S0 = x_i[2]; 
  int n_days = x_i[3]; 
  int time_seed_omicron = x_i[4];

  // define real data // 

  real sigma = x_r[1];
  real gamma = x_r[2]; 
  
  
  // define parameters need to solve model //
  
  real beta = theta[1];
  real rho  = theta[2];

  // define compartments // 
  real S;
  real E; 
  real I; 
  real Q; 
  real R;  
  
  
  // define states 
  
  real dS_dt ; 
  real dE_dt ;
  real dI_dt ;
  real dQ_dt ;
  real dR_dt ;     
  
  // define force of infection 
  real FOI;
  

// initial conitions // 
 
 S = y[1]; 
 E = y[2];
 I = y[3];
 Q = y[4]; 
 R = y[5];


  
  // calc FOI //
  
 FOI  = beta * I / S0;
  
  // ODE //
  
 dS_dt =  - FOI * S;
 dE_dt =  FOI * S - sigma * E;
 dI_dt =  (1-rho) * sigma * E - gamma * I;
 dQ_dt =  rho * sigma * E - gamma * Q ;
 dR_dt =  gamma * (I + Q); 
  
 return{dS_dt ,dE_dt,dI_dt,dQ_dt, dR_dt} ; 
  }
}
  data {
  int<lower = 1> n_ts;        // number of days observed
  int<lower = 1> n_recov;     // recovered 
  int<lower = 1> n_pop;       // population 
  int<lower = 1> n_data;      // number of data points to fit to   
  real<lower = 0> sigma;      // progression rate
  real<lower = 0> gamma;      // recovery rate
  int y[n_data];              // data, reported incidence each day 
  real ts [n_ts];             // 1:n_days 
  int <lower = 0> time_seed_omicron; // index when to fit the model to data

  }
  
  
  transformed data {

  int S0 = n_pop - n_recov;

  // any data we want to input to to our solver function 
  real x_r[2] = {sigma, gamma}; 
  int  x_i[4] = { n_recov, S0,n_ts,time_seed_omicron};

  }
  
  // this is what we want to estimate 
  
  parameters { 
  real<lower = 0> beta;
  real<lower = 0> I0;
  real<lower = 0, upper = 1> rho;
  real<lower = 0> k;

  }
  
  transformed parameters{
    
  // solution from the ODE solver 
  real y_hat[n_ts,5]; 
  
  // any parameters we want to feed to the solver 
  real theta[2] = {beta , rho}; 
  
  // initial conditions for the solver 
  real init[5]  = {S0 - I0,0, I0, 0,  n_recov };
  
  // reported incidence 
  real lambda[n_ts]; 
  real lambda_fit[n_data]; 
  
  // ODE solver 
  y_hat = integrate_ode_rk45(SIR, init, 0, ts, theta, x_r, x_i); 
  
  // calculate reported incidence for the time interval we want to fit to 
  // i.e., we allow for an initial month to seed the model 
  
 for (t in 1:n_ts)  lambda[t] = rho * sigma *  y_hat[t,2] ;

 // adding a very small number to our liklihood stops issues with the sample rejecting initial 
 // values due log probability not being >0. 
 for (t in time_seed_omicron:n_ts) lambda_fit[(t-time_seed_omicron+1)] = lambda[t] + 0.0001;
  
}
  model {
    
 // likelihood 
 
 target += neg_binomial_2_lpmf(y | lambda_fit,k);
   
 // priors
  
  
  beta ~ lognormal(0.8,0.5);
  I0   ~ normal(1,5); 
  rho  ~ beta(1,1);
  k    ~ exponential(0.01);
  }
  
  generated quantities {
  
  // basic reproduction number
  real R_0 = ((1-rho) * beta ) / gamma ; 


  }
