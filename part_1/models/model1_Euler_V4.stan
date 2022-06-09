  data { // Input data
  
  int<lower = 1> n_data;      // number of days observed
  int<lower = 1> n_ts;        // number of model time steps
  int<lower = 1> n_recov;     // recovered 
  int<lower = 1> n_pop;       // population 
 
  real<lower = 0> sigma;      // progression rate
  real<lower = 0> gamma;      // recovery rates 
  
  int y[n_data];                     // data, reported incidence each day 
  int <lower = 0> time_seed_omicron; // index when to fit the model to data 
  
  int scale_time_step ; // amount to reduce time step by when solving 
  }
  
  transformed data{
    real sigma_scale = sigma / scale_time_step;
    real gamma_scale = gamma / scale_time_step;
    int n_ts_scale = n_ts * scale_time_step ;
  }
  
  parameters { // this is what we want to estimate 
  real<lower =0> beta;
  real<lower =0> I0;
  real<lower =0, upper = 1> rho; 
  real<lower =0> k; 
  }
  
  
  transformed parameters{
    
  // compartments 
  
  real S[n_ts_scale+1];
  real E[n_ts_scale+1]; 
  real I[n_ts_scale+1]; 
  real Q[n_ts_scale+1]; 
  real R[n_ts_scale+1]; 

// Force of infection 
  real FOI[n_ts_scale];
  
// incidence 
  real lambda[n_ts_scale];

// scale beta so its a rate per scaled time step, rather than days

  real beta_scale = beta / scale_time_step;


// initial conitions 
  S[1] = n_pop -n_recov; 
  E[1] = 0;
  I[1] = I0;
  Q[1] = 0;
  R[1] = n_recov;


//// simulate ////

for(t in 1:n_ts_scale){
  
  FOI[t] = beta_scale * I[t] / (n_pop -n_recov); 
  
  S[t+1] = S[t] - FOI[t] * S[t]; 
  E[t+1] = E[t] + FOI[t] * S[t] - sigma_scale * E[t];
  I[t+1] = I[t] + (1-rho) * sigma_scale * E[t] - gamma_scale * I[t];
  Q[t+1] = Q[t] + rho * sigma_scale * E[t] - gamma_scale * Q[t]; 
  R[t+1] = R[t] + gamma_scale * (I[t] + Q[t]);
  
  lambda[t] = rho * sigma_scale * E[t]  ;
  }
  }
  
  

  model {
   
   real lambda_days[n_ts];
   real lambda_fit[n_data];
  
  // As we reduced the time step to solver our model over, 
  // we need aggregate from our reduced time step to incidence / day
  
  
  // we are going to use a for loop to sum over each day 
 
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_ts){
  ind = index + (scale_time_step-1);

  lambda_days[t] =  sum(lambda[index:ind]);

  index = index + scale_time_step;
}


  
for (t in time_seed_omicron:n_ts) lambda_fit[(t-time_seed_omicron+1)] = lambda_days[t];
  
  
    
 // likelihood 
 
 target += neg_binomial_2_lpmf(y | lambda_fit, k);
   
 // priors
  
  beta ~ lognormal(0.8,0.5);
  I0 ~ normal(1,50); 
  rho ~ beta(1,4);
  k ~ exponential(0.01);
  

  }

  
  generated quantities {
  // basic reproduction number
  real R_0 = ((1-rho) * beta ) / gamma ; 
  
   real lambda_days[n_ts];

  // As we reduced the time step to solver our model over, 
  // we need aggregate from our reduced time step to incidence / day
  
  
  // we are going to use a for loop to sum over each day 
 
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_ts){
  ind = index + (scale_time_step-1);

  lambda_days[t] =  sum(lambda[index:ind]);

  index = index + scale_time_step;
}


  }