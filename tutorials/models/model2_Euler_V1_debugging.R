# Data blocl 

n_var = 2                                # number of varaints               
n_ts = length(all_dates_mv)              # no. time steps 
n_pop = n_pop                            # population  
n_recov = n_recov_mv                     # recovered population 
y_D = multi_var_sim_inc[[1]]$rep_inc_D   # incidence data 
y_O = multi_var_sim_inc[[2]]$rep_inc_O   # incidence data
n_data_D = nrow(multi_var_sim_inc[[1]])  # no. data 
n_data_O = nrow(multi_var_sim_inc[[2]])  # no. data 
sigma = sigma                            # latent rate  
gamma = gamma                            # recovery rate 
time_seed_O =                            # index seed omicron
  which(all_dates_mv == date_seed_O)
time_fit_D = 
  which(all_dates_mv == date_fit_D)
time_fit_O = 
  which(all_dates_mv == date_fit_O)
time_int_start = which(all_dates_mv == date_int[1])
time_int_end = which(all_dates_mv == date_int[2])
scale_time_step = 5                                    # amount to reduce time step,
nu  = nu

# Parameter block 

 beta = c(2,2)
 I0 = c(1,1)
 rho = 0.2
 omega = 0.5 
 k = 0.01
 epsilon = 1/90

# Transformed data block 
 
  sigma_scale = sigma / scale_time_step;
  gamma_scale = gamma / scale_time_step;
  n_ts_scale  = n_ts * scale_time_step ;
  time_seed_O_scale = time_seed_O * scale_time_step; 
  time_int_scale_start = time_int_start * scale_time_step;
  time_int_scale_end = time_int_end * scale_time_step;
  nu_scale = nu / scale_time_step;
  
# Transformed parameters block 

S = Q = R = SO = rep(0, n_ts_scale+1)
E = I =matrix(0, nrow = n_ts_scale+1,n_var )
FOI = beta_scale = beta2 = lambda =matrix(0, nrow = n_ts_scale,n_var )


for(t in 1: n_ts_scale){
  for(i in 1:n_var)  beta2[t,i] = beta[i];
}

for(t in time_int_scale_start:time_int_scale_end){
  for(i in 1:n_var)  beta2[t,i] = beta[i] * (1-omega);
}

for(t in 1:n_ts_scale){  
  for(i in 1:n_var) beta_scale[t,i] = beta2[t,i] / scale_time_step;
}

epsilon_scale = epsilon / scale_time_step ; 

#  initial conitions 
S[1] = n_pop - n_recov; 
for(i in 1:n_var) E[1,i] = 0;
I[1,1] = I0[1];
I[1,2] = 0 ; 
Q[1] = 0;
R[1] = n_recov;
SO[1] = 0;

#  simulate
  
  for(t in 1:n_ts_scale){
    
    I[time_seed_O_scale,2] = I0[2];
    
    for(i in 1:n_var) FOI[t,i] = beta_scale[t,i] * I[t,i] / (n_pop -n_recov); 
    
    S[t+1]   =  S[t] - sum(FOI[t,]) * S[t] - nu_scale* S[t]; 
    E[t+1,1] = E[t,1] + FOI[t,1] * S[t] - sigma_scale * E[t,1];
    E[t+1,2] = E[t,2] + FOI[t,2] * (S[t] + SO[t]) - sigma_scale * E[t,2];
    
    for(i in 1:n_var) 
      I[t+1,i] = I[t,i] + (1-rho) * sigma_scale * E[t,i] - gamma_scale * I[t,i];
    Q[t+1] = Q[t] + rho * sigma_scale * sum(E[t,]) - gamma_scale * Q[t]; 
    R[t+1] = R[t] + gamma_scale * (sum(I[t,]) + Q[t]) + nu_scale *  S[t] - epsilon_scale * R[t];
    SO[t+1] = SO[t] + epsilon_scale * R[t] - FOI[t,2] * SO[t]; 
    
    for(i in 1:n_var) lambda[t,i] = rho * sigma_scale * E[t,i]  ;
  }



# Model block 

lambda_days = matrix(0, nrow = n_ts, ncol = n_var)
lambda_fit_D = rep(NA,n_data_D)
lambda_fit_O = rep(NA,n_data_O)
 

index = 1;

for (t in 1:n_ts){
  ind = index + (scale_time_step-1);
  
  for(i in 1:n_var) lambda_days[t,i] =  sum(lambda[index:ind,i]);
  
  index = index + scale_time_step;
}

#  discard 30 days before fitting Delta 
for (t in time_fit_D:(time_fit_O-1)){
  lambda_fit_D[t - time_fit_D+1] = lambda_days[t,1];
}

#  discard 30 days before fitting Omicron   
for (t in time_fit_O:n_ts) {
  lambda_fit_O[(t-time_fit_O+1)] =  lambda_days[t,2];
}