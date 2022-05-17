  data {
  int<lower = 1> n_obs;       // number of days observed
  int<lower = 1> n_var; // number of variants 
  int<lower = 1> n_recov; 
  int<lower = 1> n_pop;       // population 
 
  real<lower = 0> sigma; 
  real<lower = 0> gamma; 
  
  int y[n_obs,n_var];           // data, total number of infected individuals each day
  int scale_time_step;   // time step to estimate compartments 
  
  int time_seed_omicron; 
  int time_switch;
  int time_switch_off;
  
  }
  
  parameters { // this is what we want to estimate 
  real<lower = 0> beta[n_var];
  real<lower = 0> I0[n_var];
  real<lower = 0, upper=1> rho; 
  // real<lower =0, upper =1> omega;
  }
  
  transformed parameters{
    
    // compartments 
  real S[n_obs*scale_time_step+1];
  real E[n_obs*scale_time_step+1,n_var]; 
  real I[n_obs*scale_time_step+1,n_var]; 
  real Q[n_obs*scale_time_step+1]; 
  real R[n_obs*scale_time_step+1]; 

// Force of infection 
  real FOI[n_obs*scale_time_step,n_var];
  
// incidence 
  real lambda[n_obs*scale_time_step,n_var];
  
  real beta2[n_var];
  real omega2;



// initial conitions 
  S[1] = n_pop -n_recov - sum(I0[]); 
  for(i in 1:n_var) E[1,i] = 0;
  I[1,1] = I0[1];
  I[1,2] = I0[2];
 // I[time_seed_omicron*scale_time_step,2] = I0[2];
  Q[1] = 0;
  R[1] = n_recov;


//// simulate ////


// As we estimate beta as a daily rate, we need to scale it to our chosen time step 

 for(i in 1:n_var) beta2[i] = beta[i] /  scale_time_step;
  
 for(t in 1:n_obs*scale_time_step){
  
  
  // if  (t >= time_switch*scale_time_step && t < time_switch_off*scale_time_step) {
  //   omega2 = omega;
  // }else {
  //   omega2 = 0;
  // }

  for(i in 1:n_var) FOI[t,i] =  beta2[i] * I[t,i] / n_pop; 
  
  S[t+1] = S[t] - sum(FOI[t, ]) * S[t]; 
  for(i in 1:n_var) E[t+1,i] = E[t,i] + FOI[t,i] * S[t] - sigma * E[t,i];
  for(i in 1:n_var) I[t+1,i] = I[t,i] + (1-rho) * sigma * E[t,i] - gamma * I[t,i];
  Q[t+1] = Q[t] + rho * sigma * sum(E[t,]) - gamma * Q[t]; 
  R[t+1] = R[t] + gamma * (sum(I[t,]) + Q[t]);
  
  for(i in 1:n_var) lambda[t,i] = (1-rho) * sigma * E[t,i]  ;
  }
  
  }
  
  
  model {
    
  real lambda_days[n_obs, n_var];
  
  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_obs){
  ind = index + (scale_time_step-1);

  for(i in 1:n_var)  lambda_days[t,i] =  sum(lambda[index:ind,i]);

  index = index + scale_time_step;
}


  for(i in 1:n_var) target += poisson_lpmf(y[,i] | lambda_days[,i]);
  
  //priors
  
  beta ~ lognormal(1.5,1);
  I0 ~ normal(1,1000);  
  rho ~ beta(1,1);
  // omega ~beta(1,1);

  
  
  }
  
  generated quantities {
  real R_0[n_var];      // Basic reproduction number

  real lambda_days[n_obs,n_var];
  vector[n_var*n_obs] log_lik;
  
// used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_obs){
  ind = index + (scale_time_step-1);

  for(i in 1:n_var)  lambda_days[t,i] =  sum(lambda[index:ind,i]);

  index = index + scale_time_step;
}

  
  for(i in 1:n_obs) log_lik[i] =  poisson_lpmf(y[i,1] | lambda_days[i,1]) ; 
  for(i in n_obs+1:2*n_obs) log_lik[i] =  poisson_lpmf(y[i-n_obs,2] | lambda_days[i-n_obs,2]) ; 
  
  
  for(i in 1:n_var) R_0[i] = ( (1-rho) * beta2[i]) / gamma;

  
  }
