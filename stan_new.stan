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
  real priorSigma[2];
  real priorRange[2];
  real priorNugget[2];
  int nu; // Smoothness
  vector[N] type; // location type, urban or rural
  real max_urban; //max displacement distance constraint for urban locations
  real max_rural; //max displacement distance constraint for rural locations
  real pi;
}
parameters {
  real alpha; // intercept
  real beta; // slope
  real theta[3];// parameter
  real <lower=0, upper=2*pi> angle[N]; // random direction angle
  real <lower=0, upper=1> dist[N]; //random displacement distance for urban locations
}
transformed parameters{
  vector[N] xCoor_new; //longitudes of sampled locations
  vector[N] yCoor_new; //latitudes of sampled locations 
  matrix[N,N] dMat;//new distance matrix
  // Definitions
  real stdNugget;
  real stdSpatial;
  real range;
  cov_matrix[N] covMat;
  vector[N] eta;
  //sampling the locations
  for (i in 1:N){
    if (type[i]==1){
    xCoor_new[i]=L1[i]+max_urban*dist[i]*sin(angle[i]);
    yCoor_new[i]=L2[i]+max_urban*dist[i]*cos(angle[i]);
  } else{
    xCoor_new[i]=L1[i]+max_rural*dist[i]*sin(angle[i]);
    yCoor_new[i]=L2[i]+max_rural*dist[i]*cos(angle[i]);
  }}
  // Conversion
  stdNugget = exp(theta[1]);
  stdSpatial = exp(theta[2]);
  range = exp(theta[3]);
  //distance matrix
  for(i in 1:100){
    for(j in 1:100){
      if (i==j){
      dMat[i,j]=0;
      } else {
      dMat[i,j] = sqrt(pow(xCoor_new[i]-xCoor_new[j], 2) + pow(yCoor_new[i]-yCoor_new[j], 2));
    }
  }}
  // Make covariance matrix
  for(i in 1:100){
    for(j in 1:100){
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
  
  //calculation of random distance and random angle
 target += uniform_lpdf(dist | 0, 1);
 target += uniform_lpdf(angle | 0, 2*pi);
}
generated quantities{
  real sdNugget;
  real sdSpatial;
  real rho;
  sdNugget = stdNugget;
  sdSpatial =stdSpatial;
  rho = range;
}
