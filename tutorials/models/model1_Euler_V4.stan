  data { // Input data
  
  int<lower = 1> n_var;              // number of variants
  int<lower = 1> n_data_D;           // number of days observed (Delta)
  int<lower = 1> n_data_O;           // number of days observed (Omicron) 
  int<lower = 1> n_ts;               // number of model time steps
  int<lower = 1> n_recov;            // recovered 
  int<lower = 1> n_pop;              // population 
 
  real<lower = 0> sigma;             // progression rate
  real<lower = 0> gamma;             // recovery rates 
  
  int y_D[n_data_D];                 // data, reported incidence each day (Delta)
  int y_O[n_data_O];                 // data, reported incidence each day (Omicron)
  int <lower = 0> time_seed_O;       // index when seed Omicron
  int <lower = 0> time_fit_D ;       // index number of days to run the model for before fitting (Delta)
  int <lower = 0> time_fit_O ;       // index number of days to run the model for before fitting (Omicron)
  int <lower = 0> time_int_start;    // index when to start interventions
  int <lower = 0> time_int_end;      // index when to end interventions
  int scale_time_step ; // amount to reduce time step by when solving 
  
  real <lower = 0> v[(n_ts*scale_time_step)] ;
  }
  
  transformed data{
    real sigma_scale = sigma / scale_time_step;
    real gamma_scale = gamma / scale_time_step;
    int  n_ts_scale  = n_ts * scale_time_step ;
    int  time_seed_O_scale = time_seed_O * scale_time_step; 
    int  time_int_scale_start = time_int_start * scale_time_step;
    int  time_int_scale_end = time_int_end * scale_time_step;
  }
  
  parameters { // this is what we want to estimate 
  real<lower =0> beta[n_var];
  real<lower =1> I0[n_var];
  real<lower =0, upper = 1> rho; 
  real<lower =0, upper = 1> omega;
  real<lower =0> k; 
  }
  
  
  transformed parameters{
    
  // compartments 
  
  real  S[n_ts_scale+1];
  real  E[n_ts_scale+1,n_var]; 
  real  I[n_ts_scale+1,n_var]; 
  real  Q[n_ts_scale+1]; 
  real  R[n_ts_scale+1]; 

// Force of infection 
  real  FOI[n_ts_scale,n_var];
  
// incidence 
  real  lambda[n_ts_scale, n_var];

// scale beta so its a rate per scaled time step, rather than days

  real beta_scale[n_ts_scale, n_var];
  real beta2[n_ts_scale, n_var];
  
  for(t in 1: n_ts_scale){
    for(i in 1:n_var)  beta2[t,i] = beta[i];
  }
  
   for(t in time_int_scale_start:time_int_scale_end){
    for(i in 1:n_var)  beta2[t,i] = beta[i] * (1-omega);
  }

for(t in 1:n_ts_scale){  
  for(i in 1:n_var) beta_scale[t,i] = beta2[t,i] / scale_time_step;
}

// initial conitions 
S[1] = n_pop - n_recov; 
for(i in 1:n_var) E[1,i] = 0;
I[1,1] = I0[1];
I[1,2] = 0 ; 
Q[1] = 0;
R[1] = n_recov;

//// simulate ////

  for(t in 1:n_ts_scale){
    
    I[time_seed_O_scale,2] = I0[2];
    
    for(i in 1:n_var) FOI[t,i] = beta_scale[t,i] * I[t,i] / (n_pop -n_recov); 
    
    S[t+1] = S[t] - sum(FOI[t,]) * S[t] - v[t] * S[t]; 
    for(i in 1:n_var) E[t+1,i] = E[t,i] + FOI[t,i] * S[t] - sigma_scale * E[t,i];
    for(i in 1:n_var) I[t+1,i] = I[t,i] + (1-rho) * sigma_scale * E[t,i] - gamma_scale * I[t,i];
    Q[t+1] = Q[t] + rho * sigma_scale * sum(E[t,]) - gamma_scale * Q[t]; 
    R[t+1] = R[t] + gamma_scale * (sum(I[t,]) + Q[t]) + v[t] *  S[t];
    
    for(i in 1:n_var) lambda[t,i] = rho * sigma_scale * E[t,i]  ;
  }
  
  }
  
  

  model {
   
   real lambda_days[n_ts,n_var];
   real lambda_fit_D[n_data_D];
   real lambda_fit_O[n_data_O];
   
  // As we reduced the time step to solver our model over, 
  // we need aggregate from our reduced time step to incidence / day
  
  
  // we are going to use a for loop to sum over each day 
 
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_ts){
  ind = index + (scale_time_step-1);

  for(i in 1:n_var) lambda_days[t,i] =  sum(lambda[index:ind,i]);

  index = index + scale_time_step;
}

// discard 30 days before fitting Delta 
for (t in time_fit_D:(time_fit_O-1))
lambda_fit_D[t - time_fit_D+1] = lambda_days[t,1];

// discard 30 days before fitting Omicron   
for (t in time_fit_O:n_ts) 
lambda_fit_O[(t-time_fit_O+1)] =  lambda_days[t,2];
  
  
    
 // likelihood 
  
  target += neg_binomial_2_lpmf(y_D | lambda_fit_D, k);
  target += neg_binomial_2_lpmf(y_O | lambda_fit_O, k); 
 
 
 // priors
  
  beta ~ normal(2.2,1);
  I0 ~ normal(1,50); 
  rho ~ beta(2,8);
  k ~ exponential(0.01);
  omega ~ beta(1,1);
  

  }

  
  generated quantities {
  // basic reproduction number
  real R_0[n_var] ;

 real lambda_days[n_ts, n_var];

  // As we reduced the time step to solver our model over, 
  // we need aggregate from our reduced time step to incidence / day
  
  
  // we are going to use a for loop to sum over each day 
 
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_ts){
  ind = index + (scale_time_step-1);

  for(i in 1:n_var) lambda_days[t,i] =  sum(lambda[index:ind,i]);

  index = index + scale_time_step;
}


  for(i in 1:n_var) R_0[i] = ((1-rho) * beta[i] ) / gamma ; 
  

  }