// ---- simple_pharmacokinetic ----

functions {
  vector simplepharm(real t, vector u, real ka, real CL, real Vcent) {
    real ugut = u[1];
    real ucent = u[2];
    
    vector[2] du_dt;
    du_dt[1] = - ka * ugut;
    du_dt[2] = ka * ugut - CL / Vcent * ucent;
    
    return du_dt;
  }
}

data {
  int<lower=1> N; // number of observations
  vector[N] y; // observed concentrations (e.g., c(t))
  array[N] real timepoints; //observation timepoints
  vector[2] u0; // starting values for ugut, ucent
  real t0;
}


parameters {
  // I assume here that absorption must be positive or 0 (e.g., we do not account for
  // reverse absorption in the simple model), and that similarly clearance must be 
  // positive or 0 (no reabsorption from the periphery)
  real<lower=0> ka;
  real<lower=0> CL;
  real<lower=0> Vcent;
  real<lower=0> sigma;
}

transformed parameters {
  array[N] vector[2] u = ode_rk45(simplepharm, u0, t0, timepoints, ka, CL, Vcent);
}


model {
CL ~ lognormal(log(10), 0.25);
Vcent ~ lognormal(log(35), 0.25);
ka ~ lognormal(log(2.5), 1);
sigma ~ normal(0, 1);

for (t in 1:N) {
    y[t] ~ lognormal(log(u[t][2] / Vcent), sigma);
  }
}


generated quantities {
  array[N] real ucent_pred;
  vector[N] log_lik;
  for (t in 1:N){
    ucent_pred[t]= lognormal_rng(log(u[t][2] / Vcent), sigma);
    log_lik[t] = lognormal_lpdf(y[t] | log(u[t][2] / Vcent), sigma);
  }
}

