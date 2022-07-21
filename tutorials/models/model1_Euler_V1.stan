  data { // Input data
  
  int<lower = 1> n_data;      // number of days observed
  int<lower = 1> n_ts;        // number of model time steps
  int<lower = 1> n_recov;     // recovered 
  int<lower = 1> n_pop;       // population 
 
  real<lower = 0> sigma;      // progression rate
  real<lower = 0> gamma;      // recovery rates 
  
  int y[n_data];                     // data, reported incidence each day 
  int <lower = 0> time_seed_omicron; // index when to fit the model to data 
  real <lower = 1> I0; 
  
  }
  
  parameters { // this is what we want to estimate 
  real<lower =0> beta;
  real<lower =0, upper = 1> rho; 
  real<lower =0> k; 
  }
  
  transformed parameters{
    
  // compartments 
  
  real S[n_ts+1];
  real E[n_ts+1]; 
  real I[n_ts+1]; 
  real Q[n_ts+1]; 
  real R[n_ts+1]; 

// Force of infection 
  real FOI[n_ts];
  
// incidence 
  real lambda[n_ts];
  real lambda_fit[n_data];
  

// initial conitions 
  S[1] = n_pop -n_recov; 
  E[1] = 0;
  I[1] = I0;
  Q[1] = 0;
  R[1] = n_recov;


//// simulate ////

for(t in 1:n_ts){
  
  FOI[t] = beta * I[t] / (n_pop -n_recov); 
  
  S[t+1] = S[t] - FOI[t] * S[t]; 
  E[t+1] = E[t] + FOI[t] * S[t] - sigma * E[t];
  I[t+1] = I[t] + (1-rho) * sigma * E[t] - gamma * I[t];
  Q[t+1] = Q[t] + rho * sigma * E[t] - gamma * Q[t]; 
  R[t+1] = R[t] + gamma * (I[t] + Q[t]);
  
  lambda[t] = rho * sigma * E[t]  ;
  }

  // calculate reported incidence for the time interval we want to fit to 
  // i.e., we allow for an initial month to seed the model 
  
  
for (t in time_seed_omicron:n_ts) lambda_fit[(t-time_seed_omicron+1)] = lambda[t];
  
  
  }
  
  

  model {
    
 // likelihood 
 
 target += neg_binomial_2_lpmf(y | lambda_fit, k);
   
 // priors
  
  beta ~ normal(2.2,1);
  rho ~ beta(2,8);
  k ~ exponential(0.01);
  

  }

  
  generated quantities {
  // basic reproduction number
  real R_0 = ((1-rho) * beta ) / gamma ; 

  }