// source: https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html


functions {
   real FOMC_pdf(real x1,real x2,real threshold,real sigma,real xi,real alpha,real lambda){
       if(x1>threshold && x2>threshold){
        return(-((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(sigma^2*(1 + xi*(-threshold + x1)/sigma)*(1 + xi*(-threshold + x2)/sigma)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^2) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(alpha*sigma^2*(1 + xi*(-threshold + x1)/sigma)*(1 + xi*(-threshold + x2)/sigma)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^2));
  }
  
    if(x1>threshold && x2<=threshold){
      return(((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(sigma*(1 + xi*(-threshold + x1)/sigma)*((1/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))));
  }
    if(x1<=threshold && x2>threshold){
      return(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1 /lambda)^(-1/alpha) + ((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(sigma*(1 + xi*(-threshold + x2)/sigma)*((1/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha))));
  }
    else{    
       return(1-(((1/lambda)^(-1/alpha) + (1/lambda)^(-1/alpha))^(alpha)));
  }
  
}
//real FOMC_pdf(real x1,real x2,real threshold,real sigma,real xi,real alpha,real lambda){
//    real z1 = 1 + xi*(x1 - threshold)*(inv(sigma));
//    real z2 = 1 + xi*(x2- threshold)*(inv(sigma));
//    real inv_xi = inv(xi);
//    real inv_alpha = inv(alpha); 
//    if(x1>threshold && x2>threshold){
//      return(-(z1^inv_xi*inv(lambda))^(-inv_alpha)*(z2^inv_xi*inv(lambda))^(-inv_alpha)*((z2^inv_xi*inv(lambda))^(-inv_alpha) + (z1^inv_xi*inv(lambda))^(-inv_alpha))^alpha/(sigma^2*(z1)*(z2)*(((z2)^(1/xi)*inv(lambda))^(-inv_alpha) + ((z1)^(inv_xi)*inv(lambda))^(-inv_alpha))^2) + ((z1)^(inv_xi)*inv(lambda))^(-inv_alpha)*((z2)^(inv_xi)*inv(lambda))^(-inv_alpha)*(((z2)^(inv_xi)*inv(lambda))^(-inv_alpha) + ((z1)^(inv_xi)*inv(lambda))^(-inv_alpha))^alpha/(alpha*sigma^2*(z1)*(z2)*(((z2)^(inv_xi)*inv(lambda))^(-inv_alpha) + ((z1)^(inv_xi)*inv(lambda))^(-inv_alpha))^2));
//  }
//  
//    if(x1>threshold && x2<=threshold){
//      return(((z1)^(inv_xi)*inv(lambda))^(-inv_alpha)*(inv(lambda)^(-inv_alpha) + ((z1)^(1/xi)*inv(lambda))^(-inv_alpha))^alpha/(sigma*(z1)*(inv(lambda)^(-inv_alpha) + ((z1)^(inv_xi)*inv(lambda))^(-inv_alpha))));
//  }
//    if(x1<=threshold && x2>threshold){
//      return(((z2)^(inv_xi)*inv(lambda))^(-inv_alpha)*(inv(lambda)^(-inv_alpha) + ((z2)^(1/xi)*inv(lambda))^(-inv_alpha))^alpha/(sigma*(z2)*(inv(lambda)^(-inv_alpha) + ((z2)^(inv_xi)*inv(lambda))^(-inv_alpha))));
//  }
//    else{    
//       return(1-((inv(lambda)^(-inv_alpha) + inv(lambda)^(-inv_alpha))^(alpha)));
//  }
//  
//}
//
//


real gpd_lpdf(real y,real ymin,real sigma,real xi, real lambda){
  if(y>=ymin)
  { real inv_xi = inv(xi);
    if (xi<0 && (y-ymin)/sigma > -inv_xi)
       reject("xi<0 and max(y-ymin)/sigma > -1/xi; found xi, sigma =", xi, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (fabs(xi) > 1e-15)
      return log(lambda)-(1+inv_xi)*(log1p((y-ymin) * (xi/sigma))) -log(sigma);
    else
     return log(lambda)-(y-ymin)/sigma -log(sigma); // limit xi->0
     }
  else{
    return log((1-lambda)/ymin);
  }
}



real gpareto_lpdf(vector y, real ymin, real xi, real sigma) {
  // generalised Pareto log pdf 
  int N = rows(y);
  real inv_xi = inv(xi);
  if (xi<0 && max(y-ymin)/sigma > -inv_xi)
    reject("xi<0 and max(y-ymin)/sigma > -1/xi; found xi, sigma =", xi,",", sigma);
  if (sigma<=0)
    reject("sigma<=0; found sigma =", sigma);
  if (fabs(xi) > 1e-15)
    return -(1+inv_xi)*sum(log1p((y-ymin) * (xi/sigma))) -N*log(sigma);
  else
    return -sum(y-ymin)/sigma -N*log(sigma); // limit xi->0
}


real Hmodel_lpdf(vector x,real threshold,real scale,real xi,real lambda,real alpha){
  int N = rows(x);
  real likelihood = gpd_lpdf(x[1] | threshold,scale,xi, lambda);
  real num= 0;
  real denom =0;
  for (i in 1:(N-1)) {
    num = log(FOMC_pdf(x[i],x[i+1], threshold,scale,xi, alpha, lambda)) +num;
    denom = gpd_lpdf(x[i] | threshold,scale,xi,lambda) + denom;
  }
  likelihood = likelihood+ num -denom;
  return likelihood;
}

}
data {
  int<lower=0> S;
  int<lower=0> M;
  int<lower=0> r[M];
  matrix<lower=0>[r[1],S] y_dec;
  matrix<lower=0>[r[2],S] y_jan;
  matrix<lower=0>[r[3],S] y_feb;
  matrix<lower=0>[r[4],S] y_mar;
  matrix<lower=0>[r[5],S] y_apr;
  matrix<lower=0>[r[6],S] y_may;
  matrix<lower=0>[r[7],S] y_jun;
  matrix<lower=0>[r[8],S] y_jul;
  matrix<lower=0>[r[9],S] y_aug;
  matrix<lower=0>[r[10],S] y_sep;
  matrix<lower=0>[r[11],S] y_oct;
  matrix<lower=0>[r[12],S] y_nov;
  matrix<lower=0>[M,S] threshold;
  matrix<lower=0>[M,S] lambda;
  matrix<lower=0>[M,S] maxim;
}

  // int<lower=0> Nt;
  // vector<lower=ymin>[Nt] yt;


parameters {
  real a_sigma;
  real a_xi;
  real<lower=0.005,upper=0.020> phi_sigma;
  real<lower=0.01,upper=0.15> phi_xi;
  real<lower=0.005,upper=0.020> tau_sigma;
  real<lower=0.01,upper=0.15> tau_xi;
  vector[S] eps_sigma;
  vector[S] eps_xi;
  vector[M] gamma_sigma;
  vector[M] gamma_xi;
  vector<lower=0,upper=1>[S] alpha;
}

 // matrix[M,S] xi;

transformed parameters{
  matrix<lower=0>[M,S] sigma; 
  matrix[M,S] logsigma; 
  matrix<lower=0>[M,S] xi_star; 
  matrix[M,S] xi;
  real<lower=0> inv_phi_sigma;
  real<lower=0> inv_phi_xi;
  real<lower=0> inv_tau_sigma;
  real<lower=0> inv_tau_xi;
  inv_phi_sigma = (1/phi_sigma);
  inv_phi_xi = (1/phi_xi);
  inv_tau_sigma = (1/tau_sigma);
  inv_tau_xi = (1/tau_xi);
  
  //parameters  
  for(i in 1:M){
      for(j in 1:S){
        logsigma[i,j] = gamma_sigma[i] + eps_sigma[j];
        sigma[i,j] = exp(logsigma[i,j]);
        //xi_star[i,j] = gamma_xi[i]+ eps_xi[j] + sigma[i,j]/(maxim[i,j]-threshold[i,j]);
        //xi[i,j] = xi_star[i,j]- sigma[i,j]/(maxim[i,j]-threshold[i,j]);
        xi[i,j] = gamma_xi[i] + eps_xi[j];
        xi_star[i,j] = xi[i,j] + sigma[i,j]/(maxim[i,j]-threshold[i,j]);
      }
    }
}
  



model{
  //layer 3
  a_sigma ~ normal(0,1000); //100 if it doesn't work
  a_xi ~ normal(0,1000);
  
  phi_sigma ~ uniform(0.005,0.020);
  phi_xi ~ uniform(0.01,0.15);
  
  tau_sigma ~ uniform(0.005,0.020);
  tau_xi ~ uniform(0.01,0.15);
  
  //layer 2 

  gamma_sigma ~ normal(0,inv_tau_sigma);
  gamma_xi ~ normal(0,inv_tau_xi);
  
  
  eps_sigma ~ normal(a_sigma,inv_phi_sigma);
  eps_xi ~ normal(a_xi,inv_phi_xi);  

  //layer 1

  for(i in 1:S){
    y_dec[,i] ~ Hmodel(threshold[1,i], sigma[1,i], xi[1,i],lambda[1,i], alpha[i]);
    y_jan[,i] ~ Hmodel(threshold[2,i], sigma[2,i], xi[2,i],lambda[2,i], alpha[i]);
    y_feb[,i] ~ Hmodel(threshold[3,i], sigma[3,i], xi[3,i],lambda[3,i], alpha[i]);
    y_mar[,i] ~ Hmodel(threshold[4,i], sigma[4,i], xi[4,i],lambda[4,i], alpha[i]);
    y_apr[,i] ~ Hmodel(threshold[5,i], sigma[5,i], xi[5,i],lambda[5,i], alpha[i]);
    y_may[,i] ~ Hmodel(threshold[6,i], sigma[6,i], xi[6,i],lambda[6,i], alpha[i]);
    y_jun[,i] ~ Hmodel(threshold[7,i], sigma[7,i], xi[7,i],lambda[7,i], alpha[i]);
    y_jul[,i] ~ Hmodel(threshold[8,i], sigma[8,i], xi[8,i],lambda[8,i], alpha[i]);
    y_aug[,i] ~ Hmodel(threshold[9,i], sigma[9,i], xi[9,i],lambda[9,i], alpha[i]);
    y_sep[,i] ~ Hmodel(threshold[10,i], sigma[10,i], xi[10,i],lambda[10,i], alpha[i]);
    y_oct[,i] ~ Hmodel(threshold[11,i], sigma[11,i], xi[11,i],lambda[11,i], alpha[i]);
    y_nov[,i] ~ Hmodel(threshold[12,i], sigma[12,i], xi[12,i],lambda[12,i], alpha[i]);
  }
}