  data {
  int<lower = 1> n_obs;       // number of days observed
  int<lower = 1> n_recov; 
  int<lower = 1> n_pop;       // population 
 
  real<lower = 0> sigma; 
  real<lower = 0> gamma; 
  
  int y[n_obs];           // data, total number of infected individuals each day
  int scale_time_step;   // time step to estimate compartments 
  
  }
  
  parameters { // this is what we want to estimate 
  real<lower = 0> beta;
  real<lower = 0> I0;
  }
  
  transformed parameters{
    
    // compartments 
  real S[n_obs*scale_time_step+1];
  real E[n_obs*scale_time_step+1]; 
  real I[n_obs*scale_time_step+1]; 
  real R[n_obs*scale_time_step+1]; 

// Force of infection 
  real FOI[n_obs*scale_time_step];
  
// incidence 
  real lambda[n_obs*scale_time_step];
  
  real beta2;



// initial conitions 
  S[1] = n_pop -n_recov - I0; 
  E[1] = 0;
  I[1] = I0;
  R[1] = n_recov;


//// simulate ////


// As we estimate beta as a daily rate, we need to scale it to our chosen time step 

  beta2 = beta /  scale_time_step;
  
for(t in 1:n_obs*scale_time_step){
  
  FOI[t] = beta2 * I[t] / n_pop; 
  
  S[t+1] = S[t] - FOI[t] * S[t]; 
  E[t+1] = E[t] + FOI[t] * S[t] - sigma * E[t];
  I[t+1] = I[t] + sigma * E[t] - gamma * I[t];
  R[t+1] = R[t] + gamma * I[t];
  
  lambda[t] = sigma * E[t]  ;
  }
  
  

  }
  
  
  model {
    
  real lambda_days[n_obs];
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (i in 1:n_obs){
  ind = index + (scale_time_step-1);

    lambda_days[i] =  sum(lambda[index:ind]);

  index = index + scale_time_step;
}

target += poisson_lpmf(y | lambda_days);
  
  //priors
  
  beta ~ lognormal(1.5,1);
  I0 ~ normal(1,1000);  

  
  
  }
  
  generated quantities {
  real R_0;      // Basic reproduction number
  real gamma_days;

  
  real lambda_days[n_obs];
    vector[n_obs] log_lik;
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (i in 1:n_obs){
  ind = index + (scale_time_step-1);

    lambda_days[i] =  sum(lambda[index:ind]);

  index = index + scale_time_step;
}
  
  for(i in 1:n_obs) log_lik[i] =  poisson_lpmf(y[i] | lambda_days[i]) ; 
  
  R_0 = beta2 / gamma;

  
  }
