functions{
  real gev_lpdf(vector y, real mu, real sigma, real xi) {
    vector[rows(y)] z;
    vector[rows(y)] lp;
    int N;
    N = rows(y);
    for(n in 1:N){
      z[n] = 1 + (y[n] - mu) * xi / sigma;
      lp[n] = (1 + (1 / xi)) * log(z[n]) + pow(z[n], -1/xi);
    }
    return -sum(lp) - N * log(sigma);
  }
  
  real ggev_lpdf(real y, real mu, real sigma, real xi) {
    real z;
    real lp;
    z = 1 + (y - mu) * xi / sigma;
    lp = (1 + (1 / xi)) * log(z) + pow(z, -1/xi);
    return -lp - log(sigma);
  }
}

data
{
  int<lower = 0> N;       // number of extremes in wind speed
  vector[N] y;     // vector of extremes
  real mu_mean;
  real mu_var;
  real sigma_mean;
  real sigma_var;
  real xi_mean; 
  real xi_var;
}

transformed data {
  real min_y;
  real max_y;
  real sd_y;
  min_y = min(y);
  max_y = max(y);
  sd_y = sd(y);
}

parameters
{
  real<lower = 0> sigma; //scale
  real xi; //shape
  real mu; 		// location 

}

model
{
  
  sigma ~ normal(sigma_mean, sigma_var);
  mu ~ normal(mu_mean, mu_var); // mean and std dev
  xi ~ normal(xi_mean, xi_var);  // iid sample

  // 
  // Likelihood:
  y ~ gev(mu, sigma, xi);

}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = ggev_lpdf(y[n] | mu, sigma, xi);
  }
}
