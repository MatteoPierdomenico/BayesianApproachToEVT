// source: https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html


functions {
  real FOMC_pdf(real x1,real x2,real threshold,real sigma,real xi,real alpha,real lambda){
    real z1 = 1 + xi*(x1 - threshold)/sigma;
    real z2 = 1 + xi*(x2- threshold)/sigma;
    real inv_xi = 1/xi;
    real inv_alpha = 1/alpha; 
    real inv_lambda = 1/lambda;
    if(x1>threshold && x2>threshold){
      return (-(inv_lambda * z1^inv_xi)^(-inv_alpha)*(inv_lambda * z2^inv_xi)^(-inv_alpha)*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
              + (inv_lambda*z1^inv_xi)^(-inv_alpha))^alpha/(sigma^2*z1*z2*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
              + (inv_lambda*z1^inv_xi)^(-inv_alpha))^2) + (inv_lambda*z1^inv_xi)^(-inv_alpha)*(inv_lambda*z2^inv_xi)^(-inv_alpha)*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
              + (inv_lambda*z1^inv_xi)^(-inv_alpha))^alpha/(alpha*sigma^2*z1*z2*((inv_lambda*z2^inv_xi)^(-inv_alpha) + (inv_lambda*z1^inv_xi)^(-inv_alpha))^2));
    }
    
    if(x1>threshold && x2<=threshold){
      return ((1/(lambda*sigma))*(z1^(inv_xi-1))*((inv_lambda*z1^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z1^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)));
      
    }
    if(x1<=threshold && x2>threshold){
      return((1/(lambda*sigma))*(z2^(inv_xi-1))*((inv_lambda*z2^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z2^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)));
      
    }
    else{    
      return(1-((inv_lambda^(-1/alpha) + inv_lambda^(-1/alpha))^(alpha)));
      
    }
    
  }
  


real gpd_pdf(real x,real threshold,real scale,real xi, real lambda){
  if(x>=threshold)
  {if (fabs(xi) > 1e-15)
    return lambda*(scale^(-1)) * (1 + (xi * (x - threshold))/scale)^(-(+1/xi + 1));
    else
      return  lambda*(1/scale)*exp(-(x-threshold)/scale);}
  else{
    return (1-lambda)/threshold;
  }
}

real gpd_lpdf(real x,real threshold,real scale,real xi, real lambda){
  return log(gpd_pdf(x,threshold, scale,xi, lambda)); 
}


real gpareto_lpdf(vector y, real ymin, real xi, real sigma) {
  // generalised Pareto log pdf 
  int N = rows(y);
  real inv_xi = inv(xi);
  if (xi<0 && max(y-ymin)/sigma > -inv_xi)
    reject("xi<0 and max(y-ymin)/sigma > -1/xi; found xi, sigma =", xi, sigma);
  if (sigma<=0)
    reject("sigma<=0; found sigma =", sigma);
  if (fabs(xi) > 1e-15)
    return -(1+inv_xi)*sum(log1p((y-ymin) * (xi/sigma))) -N*log(sigma);
  else
    return -sum(y-ymin)/sigma -N*log(sigma); // limit xi->0
}


real Hmodel_lpdf(vector x,real threshold,real scale,real xi,real lambda,real alpha){
  int N = rows(x);
  real likelihood = log(gpd_pdf(x[1], threshold,scale,xi, lambda));
  real num= 0;
  real denom =0;
  for (i in 1:(N-1)) {
    num = log(FOMC_pdf(x[i],x[i+1],  threshold,scale,xi, alpha, lambda)) +num;
    denom = log(gpd_pdf(x[i],  threshold,scale,xi,lambda)) + denom;
  }
  likelihood = likelihood+ num -denom;
  return likelihood;
} 

}
data {
  int<lower=0> S;
  int<lower=0> M;
  int<lower=0> r[M];
  matrix[r[1],S] y_dec;
  matrix[r[2],S] y_jan;
  matrix[r[3],S] y_feb;
  matrix[r[4],S] y_mar;
  matrix[r[5],S] y_apr;
  matrix[r[6],S] y_may;
  matrix[r[7],S] y_jun;
  matrix[r[8],S] y_jul;
  matrix[r[9],S] y_aug;
  matrix[r[10],S] y_sep;
  matrix[r[11],S] y_oct;
  matrix[r[12],S] y_nov;
  matrix[M,S] threshold;
  matrix[M,S] lambda;
  real b;
  real c;
  real d;
  real e;
  real f;
  real g;
}
//  matrix[M,S] maxim;
  // int<lower=0> Nt;
  // vector<lower=ymin>[Nt] yt;


transformed data{
  real<lower=0> tau_a;
  tau_a =1/c;
  }


parameters {
  real a_sigma;
  real a_xi;
  real<lower=0.005,upper=0.015> phi_sigma;
  real<lower=0.05,upper=0.15> phi_xi;
  real<lower=0.005,upper=0.015> tau_sigma;
  real<lower=0.05,upper=0.15> tau_xi;
  vector[S] eps_sigma;
  vector[S] eps_xi;
  vector[M] gamma_sigma;
  vector[M] gamma_xi;
  vector<lower=0,upper=1>[S] alpha;
}

 // matrix[M,S] xi;

transformed parameters{
  matrix[M,S] sigma; 
  matrix[M,S] xi;
  real<lower=0> inv_phi_sigma;
  real<lower=0> inv_phi_xi;
  real<lower=0> inv_tau_sigma;
  real<lower=0> inv_tau_xi;
  inv_phi_sigma =(1/phi_sigma);
  inv_phi_xi =(1/phi_xi);
  inv_tau_sigma =(1/tau_sigma);
  inv_tau_xi =(1/tau_xi);
  for(i in 1:M){
      for(j in 1:S){
        sigma[i,j] = exp(gamma_sigma[i]+ eps_sigma[j]);
        xi[i,j] = gamma_xi[i]+ eps_xi[j];
      }
    }
}
  //for(i in 1:M){
  //  for(j in 1:S){
  //  xi[i,j] = xi_star[i,j] - exp(logsigma[i,j])/(maxim[i,j]-threshold[i,j]);
  //  }
  //}



model{
  //layer 3
  
  a_sigma ~ normal(b,tau_a);
  a_xi ~ normal(b,tau_a);
  
  phi_sigma ~ uniform(0.005,0.015);
  phi_xi ~ uniform(0.05,0.15);
  
  tau_sigma ~ uniform(0.005,0.015);
  tau_xi ~ uniform(0.05,0.15);
  
  //layer 2 
  for(i in 1:M){
    gamma_sigma[i] ~ normal(b,inv_tau_sigma);
    gamma_xi[i] ~ normal(b,inv_tau_xi);
  }
  for(i in 1:S){
    eps_sigma[i] ~ normal(a_sigma,inv_phi_sigma);
    eps_xi[i] ~ normal(a_xi,inv_phi_xi);
  }
  
  
  
  for(i in 1:S){ 
    alpha[i] ~ uniform(0,1);

  }
  //layer 1

  for(i in 1:S){
    target += Hmodel_lpdf(y_dec[,i]| threshold[1,i], sigma[1,i], xi[1,i],lambda[1,i], alpha[i]);
    target += Hmodel_lpdf(y_jan[,i]| threshold[2,i], sigma[2,i], xi[2,i],lambda[2,i], alpha[i]);
    target += Hmodel_lpdf(y_feb[,i]| threshold[3,i], sigma[3,i], xi[3,i],lambda[3,i], alpha[i]);
    target += Hmodel_lpdf(y_mar[,i]| threshold[4,i], sigma[4,i], xi[4,i],lambda[4,i], alpha[i]);
    target += Hmodel_lpdf(y_apr[,i]| threshold[5,i], sigma[5,i], xi[5,i],lambda[5,i], alpha[i]);
    target += Hmodel_lpdf(y_may[,i]| threshold[6,i], sigma[6,i], xi[6,i],lambda[6,i], alpha[i]);
    target += Hmodel_lpdf(y_jun[,i]| threshold[7,i], sigma[7,i], xi[7,i],lambda[7,i], alpha[i]);
    target += Hmodel_lpdf(y_jul[,i]| threshold[8,i], sigma[8,i], xi[8,i],lambda[8,i], alpha[i]);
    target += Hmodel_lpdf(y_aug[,i]| threshold[9,i], sigma[9,i], xi[9,i],lambda[9,i], alpha[i]);
    target += Hmodel_lpdf(y_sep[,i]| threshold[10,i], sigma[10,i], xi[10,i],lambda[10,i], alpha[i]);
    target += Hmodel_lpdf(y_oct[,i]| threshold[11,i], sigma[11,i], xi[11,i],lambda[11,i], alpha[i]);
    target += Hmodel_lpdf(y_nov[,i]| threshold[12,i], sigma[12,i], xi[12,i],lambda[12,i], alpha[i]);
  }
}

