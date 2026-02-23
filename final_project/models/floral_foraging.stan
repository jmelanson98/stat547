// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=0> C;                // total number of colonies
  int<lower=0> K;                // total number of traps
  int<lower=0> O;                // total number of observations (observation for a colony at a trap, in a timepoint)
  int<lower=0> CT;               // number of colonies * number of timepoints
  int CTstarts[CT];              // start index for each colony-timepoint combo
  int CTlengths[CT];             // number of obs per colony-timepoint combo
  matrix[O,2] trap_pos;              // trap coordinates
  int colony_id[O];              // colony id for each observation
  int trap_id[O];                 // trap id for each observation
  vector[O] fq;                      // floral quality for each observation
  vector[O] lq;                   // landscape quality for each observation
  int yobs[O];                   // number of individuals observed
  //real rhomax;                   // maximum foraging distance
  array[C,2] real colonycenters;      // site specific bounds on colonies
  real deltaprior;
}

parameters {
  real <lower=0> rhomax;
  real rho0; 
  real rho1; 
  real rho2; 
  real rho3; 
  real<lower=0> sigma;
  vector[K] eps;
  
  array[C] real delta_x_raw;
  array[C] real delta_y_raw;
}


transformed parameters {
  real<lower=0> sigma_sqrt = sqrt(sigma);
  vector[K] eps_scale = eps*sigma_sqrt;
  
  // transform colony locations to their specific site
  array[C] real delta_x;
  array[C] real delta_y;

  for (c in 1:C) {
    delta_x[c] = colonycenters[c,1] + delta_x_raw[c];

    delta_y[c] = colonycenters[c,2] + delta_y_raw[c];
  }
}

model {
  // set priors
  delta_x_raw ~ normal(0,deltaprior);
  delta_y_raw ~ normal(0,deltaprior);
  rhomax ~ normal(0,1);
  rho0 ~ normal(0,0.1);
  rho1 ~ normal(0,1);
  rho2 ~ normal(0,1);
  rho3 ~ normal(0,1);
  sigma ~ normal(0, 1);
  eps ~ normal(0, 1);
  
  // calculate likelihood
  for (n in 1:CT){
    int start = CTstarts[n]; // index start
    int length = CTlengths[n]; // index length
    
    // compute lambda for each trap in that landscape
    vector[length] dis = sqrt( square(delta_x[colony_id[start]] - trap_pos[start:start+length-1,1]) +
                       square(delta_y[colony_id[start]] - trap_pos[start:start+length-1,2]) );
    vector[length] rho = rhomax*inv_logit(rho0 + 
                          rho1*fq[start:start+length-1] +
                          rho2*lq[start:start+length-1] +
                          rho3*fq[start:start+length-1].*lq[start:start+length-1]);
    vector[length] lambda_it = (-(dis^2) ./ rho) + eps_scale[trap_id[start:start+length-1]];

    
    // compute multinomial probabilities and add to target likelihood
    yobs[start:start+length-1] ~ multinomial(softmax(lambda_it));
  }
}
