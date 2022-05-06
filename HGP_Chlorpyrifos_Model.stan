functions {
  
  vector GP(int N, real[] x, real rho, real alpha, vector eta, real delta) {
    vector[N] f;
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    L_K = cholesky_decompose(K);
    f = L_K * eta;
    return(f);
  }
  
  matrix GP_K(int N, real[] x, real rho, real delta) {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, 1, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    L_K = cholesky_decompose(K);
    return(L_K);
  }
  
}

data {
  
  int<lower=0> N;
  int<lower=0> N_observed;
  int<lower=0> N_missing;
  int<lower=0> N_times;
  int<lower=0> N_times_observed;
  int<lower=0> N_times_common;
  int<lower=0> N_mix;
  int<lower=0> N_rep;
  int<lower=0> N_treat;
  int<lower=0> N_species;
  int observed[N_observed];
  int missing[N_missing];
  int<lower=1,upper=N_species> species[N];
  int<lower=1,upper=N_treat> treat[N];
  int<lower=1,upper=N_mix> mix[N];
  int<lower=1,upper=N_rep> rep[N];
  int time[N];
  int time_observed[N];
  real times_scaled[N_times];
  int times_observed[N_times_observed];
  real y_obs[N_observed];
  
}
transformed data {
  
  real delta = 1e-9;
  
}

parameters {
  
  real<lower=0> sigma[N_species];
  
  real<lower=0> rho_main[N_species, N_treat];
  real<lower=0> alpha_main[N_species, N_treat];
  vector[N_times] eta_main[N_species, N_treat];
  
  real<lower=0> rho_mix[N_species, N_treat];
  real<lower=0> alpha_mix[N_species, N_treat];
  vector[N_times_observed] eta_mix[N_species, N_treat, N_mix];
  
  real<lower=0> rho_rep[N_species, N_treat];
  real<lower=0> alpha_rep[N_species, N_treat];
  vector[N_times_observed] eta_rep[N_species, N_treat, N_rep];
}

transformed parameters {
  
  matrix[N_times_observed, N_times_observed] K_mix[N_species, N_treat];
  matrix[N_times_observed, N_times_observed] K_rep[N_species, N_treat];
  
  vector[N_times] f_main[N_species, N_treat];
  vector[N_times_observed] f_mix[N_species, N_treat, N_mix];
  vector[N_times_observed] f_rep[N_species, N_treat, N_rep];
  
  vector[N] f;
  
  for (s in 1:N_species)
    for (t in 1:N_treat) {
      f_main[s,t] = GP(N_times, times_scaled, rho_main[s,t], alpha_main[s,t], eta_main[s,t], delta);
    }
  
  for (s in 1:N_species) {
    for (t in 1:N_treat) {
      K_mix[s,t] = GP_K(N_times_observed, times_scaled[times_observed], rho_mix[s,t], delta);
      for (m in 1:N_mix) {
        f_mix[s,t,m] = alpha_mix[s,t] * (K_mix[s,t] * eta_mix[s,t,m]);
      }
    }
  }
  
  for (s in 1:N_species) {
    for (t in 1:N_treat) {
      K_rep[s,t] = GP_K(N_times_observed, times_scaled[times_observed], rho_rep[s,t], delta);
      for (r in 1:N_rep) {
        f_rep[s,t,r] = alpha_rep[s,t] * (K_rep[s,t] * eta_rep[s,t,r]);
      }
    }
  }
  
  for (i in 1:N) {
    f[i] = f_main[species[i],treat[i],time[i]] + f_mix[species[i],treat[i], mix[i], time_observed[i]] + f_rep[species[i],treat[i],rep[i], time_observed[i]];
  }

}

model {
  
  for (s in 1:N_species) {
    for (t in 1:N_treat) {
      if (s <= 4) {
        rho_main[s,t] ~ inv_gamma(6.4549, 2.0289); // P[rho < 0.167] \approx 0.01, P[rho > 1] \approx 0.01
      } else {
        rho_main[s,t] ~ inv_gamma(3.64514, 0.557312); // P[rho < 0.167] \approx 0.01, P[rho > 1] \approx 0.01
      }
      alpha_main[s,t] ~ std_normal();
      eta_main[s,t] ~ std_normal();
    }
    
    for (i in 1:N_times_common) {
      eta_main[s,1,i] - eta_main[s,2,i] ~ normal(0, 0.05);
    }
  }

  for (s in 1:N_species) {
    for (t in 1:N_treat) {
      if (s <= 4) {
        rho_mix[s,t] ~ inv_gamma(6.4549, 2.0289); // P[rho < 0.167] \approx 0.01, P[rho > 1] \approx 0.01
      } else {
        rho_mix[s,t] ~ inv_gamma(3.64514, 0.557312); // P[rho < 0.167] \approx 0.01, P[rho > 1] \approx 0.01
      }
      alpha_mix[s] ~ std_normal();
      for (m in 1:N_mix) {
        eta_mix[s,t,m] ~ std_normal();
      }
    }
  }
  
  for (s in 1:N_species) {
    for (t in 1:N_treat) {
      if (s <= 4) {
        rho_rep[s,t] ~ inv_gamma(6.4549, 2.0289); // P[rho < 0.167] \approx 0.01, P[rho > 1] \approx 0.01
      } else {
        rho_rep[s,t] ~ inv_gamma(3.64514, 0.557312); // P[rho < 0.167] \approx 0.01, P[rho > 1] \approx 0.01
      }
      alpha_rep[s,t] ~ std_normal();
      for (r in 1:N_rep)
        eta_rep[s,t,r] ~ std_normal();
    }
  }
  
  for (i in 1:N_species)
    sigma[i] ~ student_t(3, 0, 5);
  
  y_obs ~ normal(f[observed], sigma[species[observed]]);
  
}

generated quantities {
  
  real y_predict[N] = normal_rng(f, sigma[species]);
  vector[N_missing] y_missing = f[missing];
  
}
