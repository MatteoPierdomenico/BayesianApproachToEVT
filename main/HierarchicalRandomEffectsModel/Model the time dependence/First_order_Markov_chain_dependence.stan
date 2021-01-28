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

//   real FOMC_pdf(real x1,real x2,real threshold,real sigma,real xi,real alpha,real lambda){
//    real z1 = 1 + xi*(x1 - threshold)/sigma;
//    real z2 = 1 + xi*(x2- threshold)/sigma;
//    real inv_xi = 1/xi;
//    real inv_alpha = 1/alpha; 
//    //real inv_lambda = 1/lambda;
//    if(x1>threshold && x2>threshold){
//     //return (-(inv_lambda * z1^inv_xi)^(-inv_alpha)*(inv_lambda * z2^inv_xi)^(-inv_alpha)*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
//     //       + (inv_lambda*z1^inv_xi)^(-inv_alpha))^alpha/(sigma^2*z1*z2*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
//     //       + (inv_lambda*z1^inv_xi)^(-inv_alpha))^2) + (inv_lambda*z1^inv_xi)^(-inv_alpha)*(inv_lambda*z2^inv_xi)^(-inv_alpha)*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
//     //       + (inv_lambda*z1^inv_xi)^(-inv_alpha))^alpha/(alpha*sigma^2*z1*z2*((inv_lambda*z2^inv_xi)^(-inv_alpha) + (inv_lambda*z1^inv_xi)^(-inv_alpha))^2));
//    //
//    return(-(z1^inv_xi/lambda)^(-inv_alpha)*(z2^inv_xi/lambda)^(-inv_alpha)*((z2^inv_xi/lambda)^(-inv_alpha) + (z1^inv_xi/lambda)^(-inv_alpha))^alpha/(sigma^2*(z1)*(z2)*(((z2)^(1/xi)/lambda)^(-inv_alpha) + ((z1)^(inv_xi)/lambda)^(-inv_alpha))^2) + ((z1)^(inv_xi)/lambda)^(-inv_alpha)*((z2)^(inv_xi)/lambda)^(-inv_alpha)*(((z2)^(inv_xi)/lambda)^(-inv_alpha) + ((z1)^(inv_xi)/lambda)^(-inv_alpha))^alpha/(alpha*sigma^2*(z1)*(z2)*(((z2)^(inv_xi)/lambda)^(-inv_alpha) + ((z1)^(inv_xi)/lambda)^(-inv_alpha))^2));
//  }
//  
//    if(x1>threshold && x2<=threshold){
//     //return ((1/(lambda*sigma))*(z1^(inv_xi-1))*((inv_lambda*z1^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z1^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)));
//    return(((z1)^(inv_xi)/lambda)^(-inv_alpha)*((1/lambda)^(-inv_alpha) + ((z1)^(1/xi)/lambda)^(-inv_alpha))^alpha/(sigma*(z1)*((1/lambda)^(-inv_alpha) + ((z1)^(inv_xi)/lambda)^(-inv_alpha))));
//  }
//    if(x1<=threshold && x2>threshold){
//     //return((1/(lambda*sigma))*(z2^(inv_xi-1))*((inv_lambda*z2^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z2^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)));
//return(((z2)^(inv_xi)/lambda)^(-inv_alpha)*((1/lambda)^(-inv_alpha) + ((z2)^(1/xi)/lambda)^(-inv_alpha))^alpha/(sigma*(z2)*((1/lambda)^(-inv_alpha) + ((z2)^(inv_xi)/lambda)^(-inv_alpha))));
//  }
//    else{    
//     //return(1-((inv_lambda^(-1/alpha) + inv_lambda^(-1/alpha))^(alpha)));
//       return(1-(((1/lambda)^(-inv_alpha) + (1/lambda)^(-inv_alpha))^(alpha)));
//  }
//  
//}
//
//   real FOMC_pdf(real x1,real x2,real threshold,real sigma,real xi,real alpha,real lambda){
//    real z1 = 1 + xi*(x1 - threshold)/sigma;
//    real z2 = 1 + xi*(x2- threshold)/sigma;
//    real inv_xi = 1/xi;
//    real inv_alpha = 1/alpha; 
//    real inv_lambda = 1/lambda;
//    real l1 = inv_lambda*z1^inv_xi;
//    real l2 = inv_lambda*z2^inv_xi;
//    if(x1>threshold && x2>threshold){
//    return(-l1^(-inv_alpha)*l2^(-inv_alpha)*(l2^(-inv_alpha) + l1^(-inv_alpha))^alpha/(sigma^2*z1*z2*(l2^(-inv_alpha) 
//           +l1^(-inv_alpha))^2) + l1^(-inv_alpha)*l2^(-inv_alpha)*(l2^(-inv_alpha) 
//           +l1^(-inv_alpha))^alpha/(alpha*sigma^2*z1*z2*(l2^(-inv_alpha) + l1^(-inv_alpha))^2));
//  }
//  
//    if(x1>threshold && x2<=threshold){
//     //return ((1/(lambda*sigma))*(z1^(inv_xi-1))*((inv_lambda*z1^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z1^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)));
//    return(l1^(-1/alpha)*(inv_lambda^(-inv_alpha) + l1^(-inv_alpha))^alpha/(sigma*z1*(inv_lambda^(-inv_alpha) + l1^(-inv_alpha))));
//  }
//    if(x1<=threshold && x2>threshold){
//     //return((1/(lambda*sigma))*(z2^(inv_xi-1))*((inv_lambda*z2^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z2^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)));
//    return(l2^(-inv_alpha)*((inv_lambda)^(-inv_alpha) + l2^(-inv_alpha))^alpha/(sigma*z2*(inv_lambda^(-inv_alpha) + l2^(-inv_alpha))));
//  }
//    else{    
//     //return(1-((inv_lambda^(-1/alpha) + inv_lambda^(-1/alpha))^(alpha)));
//       return(1-((inv_lambda^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha)));
//  }
//  
//}


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
  int<lower=0> N;
  vector[N] y;
  real<lower=0> threshold;
  real<lower=0> lambda;
  real b;
  real c;
  real f;
  real g;
  // int<lower=0> Nt;
  // vector<lower=ymin>[Nt] yt;
}

transformed data{
  real<lower=0> tau_a;
  tau_a =1/c;
  
}

parameters {
  real a_sigma;
  real a_xi;
  real<lower=0.005,upper=0.020> phi_sigma;
  real<lower=0.01,upper=0.10> phi_xi;
  real<lower=0,upper=1> alpha;
  real<lower=0> sigma; 
  real<lower=-sigma/(max(y)-threshold)> xi;

}
transformed parameters{
  real<lower=0> tau_sigma;
  real<lower=0> tau_xi;
  tau_sigma =(1/phi_sigma);
  tau_xi =(1/phi_xi);
}

model {
  
  //layer 3
  a_sigma ~ normal(b,tau_a);
  a_xi ~ normal(b,tau_a);
  //phi_sigma ~ gamma(f,g);    
  //phi_xi ~ gamma(f,g);
  phi_sigma ~ uniform(0.005,0.020);
  phi_xi ~ uniform(0.01,0.10);
  //layer 2 
  
  sigma ~ normal( a_sigma,tau_sigma);
  
  //sigma ~ normal(0,100);
  xi ~  normal( a_xi,tau_xi);
    
  //xi ~ normal(0,10);

  //sigma ~ lognormal(a_sigma,1/phi_sigma);
  //xi ~ normal(a_xi,1/phi_xi);
  alpha ~ uniform(0,1);
  //layer 1
  
  //for(i in 1:M){
  //  target += Hmodel_lpdf(y[,i]| threshold[i], sigma[i], xi[i],lambda[i], alpha[i]);
  //}
  
  //for(i in 1:M){
  //  target += gpareto_lpdf(y[,i]| threshold[i], xi[i], sigma[i]);
  //}
  //
  //y ~ Hmodel(threshold,sigma,xi, lambda, alpha);
  y ~ Hmodel(threshold,sigma,xi, lambda, alpha);
  //
   
}