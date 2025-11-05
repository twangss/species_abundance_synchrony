data {
  int<lower=0> N; // number of years
  int<lower=0> M; // number of time series
  array [M]int<lower=0> states;
  vector[M] add_up; // vector to add up for total abundance trend
  // int<lower=0> states[M]; // vector assigning time series to states
  int<lower=0> S; // number of states
  array [M]int<lower=0> obsVariances;  // observation variance map
  // int<lower=0> obsVariances[M];  // observation variance map
  int<lower=0> n_obsvar;
  array [S+1]int<lower=0> proVariances;  // process variance map
  // int<lower=0> proVariances[S+1];  // process variance map
  int<lower=0> n_provar;
  array [S+1]int<lower=0> trends;  // trend map
  // int<lower=0> trends[S+1];  // trend map
  int<lower=0> n_trends;
  int<lower=0> n_pos; // number of non-NA values
  array [n_pos]int<lower=0> col_indx_pos;
  array [n_pos]int<lower=0> row_indx_pos;
  int<lower=0> est_trend;
  // int<lower=0> est_sigma_process_prior;
  // vector[n_provar] sigma_process_prior; // sigma process priors
  //int<lower=0> est_B;
  int<lower=0> n_A;
  int<lower=0> est_nu;
  array [n_A+2]int<lower=0> est_A;
  vector[n_pos] y; // data
  // int y_int[n_pos];
  int family; // 1 = normal, 2 = binomial, 3 = poisson, 4 = gamma, 5 = lognormal
  
  vector[2] sigma_obs_priors; // priors 
  vector[2] sigma_proc_priors; // priors 
}
parameters {
  vector[S] x0; // initial states
  vector<lower=2>[est_nu] nu; // nu, constrainted to be > 2
  matrix <lower=-3,upper=3>[N-1,S] pro_dev;
  // vector<lower=-3,upper=3>[S] pro_dev[N-1];
  // vector[n_trends * est_trend] U;
  vector[n_trends] U;
  //matrix[S*est_B,S*est_B] B;
  vector[n_A] A; // offsets
  array [n_provar]real<lower=0> sigma_process; //real?
  array [n_obsvar]real<lower=0> sigma_obs; //real?
  // vector<lower=0,upper=1>[S] gomp_dd; // initial states
}
transformed parameters {
  real<lower=0> sigma_process_real = sigma_process[1];
  matrix[N, M] pred;
  // vector[M] pred[N];
  matrix[N,S] x; // elements accessed [N,K]
  // vector[S] x[N]; // elements accessed [N,K]
  //matrix[N, S] x;
  //matrix[S,S] Bmat;
  vector[S] Uvec;
  vector [M] Avec;
  
  for(i in 1:M) Avec[i] = 0;
  for(i in 1:n_A) Avec[est_A[i]] = A[i];
  
  for(i in 1:S) {
    if(est_trend) {
      Uvec[i] = U[trends[i]]; // map shared trends
    } else {
     Uvec[i] = 0;
    }
  }
  // for(i in 1:S) {
  //   for(j in 1:S) {
  //     if(i==j) {
  //       Bmat[i,j] = 1;
  //     } else {
  //       Bmat[i,j] = 0;
  //     }
  //   }
  // }
  //if(est_B) Bmat = B;
  
  for(s in 1:S) {x[1,s] = x0[s];}
  for(t in 2:N) {
    for(s in 1:S) {
    //x[t,] = Bmat * x[t-1,] + pro_dev[t-1,];
    //x[t,] = x[t-1,] + pro_dev[t-1,];
    x[t,s] = x[t-1,s] + pro_dev[t-1,s] * sigma_process[proVariances[s]]; // * sigma_obs[1]; 
    if(est_trend == 1) {
     x[t,s] = x[t,s] + Uvec[s];
    }
    }
  }

  // map predicted states to time series
  for(m in 1:M) {
    for(t in 1:N) {
      pred[t,m] = x[t,states[m]] + Avec[m];
    }
  }
}
model {
  //x0 ~ normal(0, 3); // initial states
  //A ~ normal(0, 3); // A offsets
  
  for(i in 1:n_obsvar) {
    sigma_obs[i] ~ inv_gamma(sigma_obs_priors[1],sigma_obs_priors[2]) ;//student_t(5,0,0.1); // observation var sigma
  }
  for(s in 1:n_provar) {
    sigma_process[s] ~ inv_gamma(sigma_proc_priors[1],sigma_proc_priors[2])  ; // inv_gamma(0.01,0.03) T[1e-4,1e2]; //inv_gamma(0.01,0.03) T[1e-3,1]; //student_t(5,0,0.1); // process var sigma
    // sigma_process[s] ~ normal(0,0.2) T[1e-3,1e2]; turned off for beta prior
    //if(est_sigma_process_prior == 1){
    //  sigma_process[s] ~ normal(sigma_process_prior[s],0.2);
    // }
  }
  
  // for(s in 1:S) {
  //    gomp_dd[s] ~ beta(1,1);//student_t(5,0,0.1); // observation var sigma
  // }
  
  //if(est_trend==1){
    for(i in 1:n_trends) {
      U[i] ~ normal(0,0.1); // optional trends
    }
  //}

  
  if(est_nu ==1) {
    nu[1] ~ gamma(2, 0.1);
    // for(s in 1:S) { // drawn from student-t distribution
    //   pro_dev[s] ~ student_t(nu[1], 0, 1); // process deviations
    // }
  for(t in 1:(N-1)) {
      for(s in 1:S) {
        pro_dev[t,s] ~ student_t(nu[1], 0, 1);
      }
    }
  } else {
    // for(s in 1:S) { // drawn from normal distribution
    //   pro_dev[s] ~ normal(0, 1); // process deviations
    // }
    for(t in 1:(N-1)) {
      for(s in 1:S) {
        pro_dev[t,s] ~ normal(0,1);
      }
    }
  }
  

  // if(est_B ==1) {
  //   for(i in 1:S) {
  //     for(j in 1:S) {
  //       if(i==j) {
  //         B[i,j] ~ uniform(0,1);
  //       } else {
  //         B[i,j] ~ normal(0,1);
  //       }
  //     }
  //   }
  // }

  // likelihood
  if(family == 1) {
    for(i in 1:n_pos) y[i] ~ normal(pred[col_indx_pos[i], row_indx_pos[i]], sigma_obs[obsVariances[row_indx_pos[i]]]);
  }
  //if(family == 2) {
    //for(i in 1:n_pos) y_int[i] ~ bernoulli_logit(pred[col_indx_pos[i], row_indx_pos[i]]);
  //}
  // if(family == 3) {
  //   for(i in 1:n_pos) y_int[i] ~ poisson_log(pred[col_indx_pos[i], row_indx_pos[i]]);
  // }
  if(family == 4) {
    for(i in 1:n_pos) y[i] ~ gamma(sigma_obs[obsVariances[row_indx_pos[i]]], sigma_obs[obsVariances[row_indx_pos[i]]] ./ pred[col_indx_pos[i], row_indx_pos[i]]);
  }
  if(family == 5) {
    for(i in 1:n_pos) y[i] ~ lognormal(pred[col_indx_pos[i], row_indx_pos[i]], sigma_obs[obsVariances[row_indx_pos[i]]]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
//   // regresssion example in loo() package
  if(family==1) for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[col_indx_pos[n], row_indx_pos[n]], sigma_obs[obsVariances[row_indx_pos[n]]]);
//   // if(family==2) for (n in 1:n_pos) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[col_indx_pos[n], row_indx_pos[n]]));
//   // if(family==3) for (n in 1:n_pos) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[col_indx_pos[n], row_indx_pos[n]]));
//   if(family==4) for (n in 1:n_pos) log_lik[n] = gamma_lpdf(y[n] | sigma_obs[obsVariances[row_indx_pos[n]]], sigma_obs[obsVariances[row_indx_pos[n]]] ./ exp(pred[col_indx_pos[n], row_indx_pos[n]]));
//   if(family==5) for (n in 1:n_pos) log_lik[n] = lognormal_lpdf(y[n] | pred[col_indx_pos[n], row_indx_pos[n]], sigma_obs[obsVariances[row_indx_pos[n]]]);
vector[N] tot_sum_x; 
for (n in 1:N) tot_sum_x[n] = dot_product(exp(pred[n,]),add_up) ; 

  // for(m in 1:M) {
  //   for(t in 1:N) {
  //     pred[t,m]*add_up[m]
  //   }
  // }
  
  
}
