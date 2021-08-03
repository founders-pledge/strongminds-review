
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] estimated_control_effect;
  vector[N] mean_change;
  vector[N] pre_n;
  vector[N] pct_attrition;
}

transformed data {
  vector[N] effect;
  effect = mean_change - estimated_control_effect;
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real theta[N];
  real<lower=0> sigma;
  real<lower=0> tau;
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // random effects with intervention type effects 
  //mean_change ~ normal(mean_change_hat, sigma);
  
  // modeling the effect sizes with a conservative prior
  //mu ~ gamma(1, 2); // what should the prior actually be?
  effect ~ normal(theta, sigma);
  theta ~ normal(mu, tau);
  mu ~ normal(5, 1); // default informative prior here for effect size?
  tau ~ cauchy(0,5);

}

