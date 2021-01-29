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
  return num;
}
}

data {
  int<lower=0> N;
  int<lower=0> M;
  matrix[N,M] y;
  vector<lower=0>[M] threshold;
  vector<lower=0>[M] lambda;
  real b;
  real c;
}

transformed data{
  real<lower=0> tau_a;
  tau_a =1/c;
  
}


parameters {
  real a_sigma;
  real a_xi;
  real<lower=0.005,upper=0.015> phi_sigma;
  real<lower=0.05,upper=0.15> phi_xi;
  vector<lower=0,upper=1>[M] alpha;
  vector[M] logsigma;
  vector[M] xi; 

}
transformed parameters{
  vector<lower=0>[M] xi_star;
  vector<lower=0>[M] sigma; 
  real<lower=0> tau_sigma;
  real<lower=0> tau_xi;
  tau_sigma =(1/phi_sigma);
  tau_xi =(1/phi_xi);
  for(i in 1:M){
    sigma[i] = exp(logsigma[i]);
    xi_star[i] = xi[i] + sigma[i]/(max(y[,i])-threshold[i]);
  }
}

  

model {
  //layer 3
  a_sigma ~ normal(b,tau_a);
  a_xi ~ normal(b,tau_a);

  phi_sigma ~ uniform(0.005,0.015);
  phi_xi ~ uniform(0.05,0.15);
  
  //layer 2 
  
  logsigma ~ normal(a_sigma,tau_sigma);
  xi ~ normal(a_xi,tau_xi);
  alpha ~ uniform(0,1);
  
  //layer 1
  
  for(i in 1:M){
    y[,i] ~ Hmodel(threshold[i], sigma[i], xi[i],lambda[i], alpha[i]);
  }
  //for(i in 1:M){
  //  for(j in 1:N){
  //    target += gpd_lpdf(y[j,i]| threshold[i], sigma[i], xi[i], lambda[i]);
  //}
  //}
}