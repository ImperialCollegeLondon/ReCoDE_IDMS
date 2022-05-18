  data {
  int<lower = 1> n_days;       // number of days observed
  int<lower = 1> n_var; // number of variants 
  int<lower = 1> n_recov; 
  int<lower = 1> n_pop;       // population 
  int<lower = 1> n_dataD; 
  int<lower = 1> n_dataO;        
 
  real<lower = 0> sigma; 
  real<lower = 0> gamma; 
  
  int y_D[n_dataD];           // data, total number of infected individuals each day
  int y_O[n_dataO];           // data, total number of infected individuals each day

  int scale_time_step;   // time step to estimate compartments 
  
  int n_ts ; // number of estimation time steps 
  
  int time_seed_omicron; 
  int time_switch;
  int time_switch_off;
  
  }
  
  parameters { // this is what we want to estimate 
  real<lower = 0> beta[n_var];
  real<lower = 0> I0[n_var];
  real<lower = 0, upper=1> rho; 
  real<lower =0, upper =1> omega;
  }
  
  transformed parameters{
    
    // compartments 
  real S[n_ts+1];
  real E[n_ts+1,n_var]; 
  real I[n_ts+1,n_var]; 
  real Q[n_ts+1]; 
  real R[n_ts+1]; 

// Force of infection 
  real FOI[n_ts,n_var];
  
// incidence 
  real lambda[n_ts,n_var];
  
  real beta2[n_ts,n_var];




// initial conitions 
  S[1] = n_pop -n_recov - I0[1]; 
  for(i in 1:n_var) E[1,i] = 0;
  for(i in 1:n_var) I[1,i] = I0[i];
  I[1,1] = I0[1];
  I[1,2] = 0;

  
  Q[1] = 0;
  R[1] = n_recov;


// As we estimate beta as a daily rate, we need to scale it to our chosen time step 
for (t in 1:n_ts){
   for(i in 1:n_var) beta2[t,i] = (beta[i] /  scale_time_step); 
}

for (t in time_switch:time_switch_off){
   for(i in 1:n_var) beta2[t,i] = (beta[i] /  scale_time_step) * (1-omega); 
}






//// simulate ////

 for(t in 1:n_ts){
  
  I[time_seed_omicron,2] = I0[2];

  for(i in 1:n_var) FOI[t,i] =  beta2[t,i] * I[t,i] / n_pop; 
  
  S[t+1] = S[t] - sum(FOI[t, ]) * S[t]; 
  for(i in 1:n_var) E[t+1,i] = E[t,i] + FOI[t,i] * S[t] - sigma * E[t,i];
  for(i in 1:n_var) I[t+1,i] = I[t,i] + (1-rho) * sigma * E[t,i] - gamma * I[t,i];
 
  Q[t+1] = Q[t] + rho * sigma * sum(E[t,]) - gamma * Q[t]; 
  R[t+1] = R[t] + gamma * (sum(I[t,]) + Q[t]);
  
  for(i in 1:n_var) lambda[t,i] = (1-rho) * sigma * E[t,i]  ;
  }
  
  }
  
  
  model {
    
  real lambda_days[n_days, n_var];
  real lambdaD [n_dataD];  
  real lambdaO [n_dataO];  

  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_days){
  ind = index + (scale_time_step-1);

  for(i in 1:n_var)  lambda_days[t,i] =  sum(lambda[index:ind,i]);

  index = index + scale_time_step;
}



 for (i in 1:n_dataD)    lambdaD[i] = lambda_days[i+(n_days-n_dataD),1] + 0.0000000001;
 for (i in 1:n_dataO)    lambdaO[i] = lambda_days[i+(n_days-n_dataO),2] + 0.0000000001;


 target += poisson_lpmf(y_D | lambdaD);
 target += poisson_lpmf(y_O | lambdaO);

  
  //priors
  
  beta ~ lognormal(1.5,1);
  I0 ~ normal(1,100);  
  rho ~ beta(1,1);
  omega ~ beta(1,1);

  
  
  }
  
  generated quantities {
  real R_0[n_var];      // Basic reproduction number

  real lambda_days[n_days, n_var];
  real lambdaD [n_dataD];  
  real lambdaO [n_dataO]; 
  
  
  vector[n_dataD+n_dataO] log_lik;

  // used for for loop
  int index;
  int ind ;

   index = 1;

  for (t in 1:n_days){
  ind = index + (scale_time_step-1);

  for(i in 1:n_var)  lambda_days[t,i] =  sum(lambda[index:ind,i]);

  index = index + scale_time_step;
}



 for (i in 1:n_dataD)    lambdaD[i] = lambda_days[i+(n_days-n_dataD),1] ;
 for (i in 1:n_dataO)    lambdaO[i] = lambda_days[i+(n_days-n_dataO),2] ;

  

  for(i in 1:n_dataD) log_lik[i] =   poisson_lpmf(y_D[i] | lambdaD[i]) ; 
  for(i in n_dataD+1:n_dataD+n_dataO) log_lik[i] =  poisson_lpmf(y_O[i-n_dataD] | lambdaO[i-n_dataD]) ; 
  
  
  for(i in 1:n_var) R_0[i] = ( (1-rho) * beta[i]) / gamma;

  
  }
