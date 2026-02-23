// DISTANCES IN KM, NOT METERS!!!
  
data {
  int<lower=1> C;                // # colonies
  int<lower=1> K;                // # traps
  int<lower=1> O;                // # observations
  int starts[C];                 // start indices for each colony
  int lengths[C];                // lengths for each colony
  matrix[O, 2] trap_pos;         // trap coordinates
  vector[O] sample_effort;      // sampling effort for each obs
  int colony_id[O];              // colony id for each obs
  int trap_id[O];                // trap id for each obs
  int y_obs[O];                  // concatenated counts for all colonies
  vector[O] yn;                 // are there any bees from that colony at that trap? 0/1
  array[C,4] real colonybounds;      // site specific bounds on colonies
  real rho;
}


parameters {

  // using hard bounds here so that init doesn't fail with - inf
  array[C] real<lower=0, upper=1> delta_x_raw;
  array[C] real<lower=0, upper=1> delta_y_raw;
}

transformed parameters {
    // transform colony locations to their specific site
  array[C] real delta_x;
  array[C] real delta_y;

  for (c in 1:C) {
    delta_x[c] = colonybounds[c,1] +
                  delta_x_raw[c] * (colonybounds[c,2] - colonybounds[c,1]);

    delta_y[c] = colonybounds[c,3] +
                  delta_y_raw[c] * (colonybounds[c,4] - colonybounds[c,3]);
  }
  
}

model {
  
  {
  
  // calculate likelihood
  for (c in 1:C) {
    
    // get some indices
    int start = starts[c]; // index start
    int length = lengths[c]; // index length
    
    // put prior on distance to traps where we find bees
    vector[length] dis = sqrt( square(delta_x[colony_id[start]] - trap_pos[start:start+length-1,1]) +
                       square(delta_y[colony_id[start]] - trap_pos[start:start+length-1,2]) );
    
    // compute visitation intensity lambda for each trap in that landscape
    vector[length] lambda_ik = -0.5*(dis / rho)^2 +
                      log(sample_effort[start:start+length-1]);
    
    // add to target likelihood
    y_obs[start:start+length-1] ~ multinomial(softmax(lambda_ik));
  }
}
}

generated quantities {
  vector[C] loglik;
    {
  
  // calculate likelihood
  for (c in 1:C) {
    
    // get some indices
    int start = starts[c]; // index start
    int length = lengths[c]; // index length
    
    // put prior on distance to traps where we find bees
    vector[length] dis = sqrt( square(delta_x[colony_id[start]] - trap_pos[start:start+length-1,1]) +
                       square(delta_y[colony_id[start]] - trap_pos[start:start+length-1,2]) );

    // compute visitation intensity lambda for each trap in that landscape
    vector[length] lambda_ik = -0.5*(dis / rho)^2 +
                      log(sample_effort[start:start+length-1]);
    
    // generate ypred values
    loglik[c] = multinomial_lpmf(y_obs[start:start+length-1] | softmax(lambda_ik));
  }
}
  
}
