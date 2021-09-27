
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  real adjusted_change[N];
  real<lower=0> change_sd[N];
  real prior_mean;
  real<lower=0> prior_sd;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu; // mean ("true") treatment effect
  real theta_tilde[N]; // per-trial treatment effect
  real<lower=0> tau; // sd of treatment effects
}

/// divergent transitions means we have to reparameterize the model

transformed parameters {
  real theta[N];
  for (n in 1:N)
    theta[n] = mu + tau * theta_tilde[n];
}


// basically what's happening in this reparameterization is this
// we say our effect size is distributed normally with mean theta and sd change_sd.
// theta  is parametrized in terms of a mean treatment effect mu, with normal errors
// those normal errors are given by theta_tilda (with a normal(0,1) prior)
// and their magnitude is given by tau, for which we have a half-cauchy prior

model {
  // priors
  mu ~ normal(prior_mean, prior_sd); // default informative prior here for effect size?
  tau ~ cauchy(0,5);
  
  // modeling the effect sizes with a conservative prior
  adjusted_change ~ normal(theta, change_sd);
  theta_tilde ~ normal(0, 1);
}

