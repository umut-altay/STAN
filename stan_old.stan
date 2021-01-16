functions {
  // PC prior defined for log-precision on log-scale
    real pc_logprec_lpdf(real par, real U, real alpha){
    real lambda = -log(alpha)/U;
    return(log(lambda/2) - par/2 - lambda*exp(-par/2));
  }
  // PC prior defined for log-precision on log-scale
  real pc_logstd_lpdf(real par, real U, real alpha){
    real lambda = -log(alpha)/U;
    return(log(lambda) + par - lambda*exp(par));
  }
  // PC prior for range
  real pc_range_lpdf(real par, real U, real alpha){
    real lambda = -log(alpha)*U;
    return(log(lambda)-par-lambda*exp(-par));
  }
}
data {
  int<lower=1> N; //the number of observations
  vector[N] x; // covariate
  vector[N] y; // response
  real L1[N]; //longitude
  real L2[N];//latitude
  int nu; // Smoothness
  real priorSigma[2];
  real priorRange[2];
  real priorNugget[2];
}

transformed data{
  matrix[N,N] dMat;
  for(i in 1:N){
    for(j in 1:N){
      dMat[i,j] = sqrt(pow(L1[i]-L1[j], 2) + pow(L2[i]-L2[j], 2));
    }
  }
}
parameters {
  real alpha; // intercept
  real beta; // slope
  real<lower=0> theta[3];// parameter
}
transformed parameters{
  // Definitions
  real stdNugget;
  real stdSpatial;
  real range;
  cov_matrix[N] covMat;
  vector[N] eta;
  // Conversion
  stdNugget = exp(theta[1]);
  stdSpatial = exp(theta[2]);
  range = exp(theta[3]);
  // Make covariance matrix
  for(i in 1:N){
    for(j in 1:N){
      if(i == j){
        covMat[i,j] = pow(stdSpatial, 2) + pow(stdNugget, 2);
      } else{
      covMat[i,j] = pow(stdSpatial,2)*modified_bessel_second_kind(nu, dMat[i,j]*(sqrt(8*nu)/range))*pow(sqrt(8*nu)/range*dMat[i,j], nu)*pow(2,1-nu)/exp(lgamma(nu)); 
      }
    }
  }
  covMat = 0.5*(covMat + covMat');
  // Linear predictor
     eta = alpha + beta*x;
 }
model{
  // Priors
  target += pc_logstd_lpdf(theta[1]|priorNugget[1], priorNugget[2]); // prior for the log-precision parameter
  target += pc_logstd_lpdf(theta[2]|priorSigma[1], priorSigma[2]);
  target += pc_range_lpdf(theta[3]|priorRange[1], priorRange[2]);
  // Model
  target +=normal_lpdf(alpha|0, 100); //prior for the intercept
  target +=normal_lpdf(beta|0, 100); //prior for the slope
  // Likelihood
 target += multi_normal_lpdf(y | eta, covMat); //was previously commented out
}
generated quantities{
  real sdNugget;
  real sdSpatial;
  real rho;
  sdNugget = stdNugget;
  sdSpatial =stdSpatial;
  rho = range;
}







