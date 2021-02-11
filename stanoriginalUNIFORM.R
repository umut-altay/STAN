## Preliminaries
# Clear the environment
rm(list = ls())

# Load libraries
library(rstan)
library(coda)
library(INLA)
library(rgdal)
#install.packages("spatialEco")
library(spatialEco)

inla.setOption(pardiso.license = "simulation_uniform/pardiso.lic")

# Set seed
set.seed(2015104)

## Functions
# covariance function
covFun = function(dMat, range, stdDev){
  Sig = inla.matern.cov(nu = 1,
                        kappa = sqrt(8*1)/range,
                        x = dMat,
                        corr = TRUE)
  Sig = stdDev^2*Sig
}

## Extract spatial design based on Kenya dataset
# Read Kenya geography
kenya.geo <-readOGR('simulation_uniform/geodata/KEGE71FL.shp')

# Remove clusters with missing locations
kenya.geo = subset(kenya.geo, kenya.geo$LONGNUM != 0)

# Prediction grid
xx = seq(min(kenya.geo$LONGNUM), max(kenya.geo$LONGNUM), length.out = 50)
yy = seq(min(kenya.geo$LATNUM), 5.5, length.out = 50) #reason for 5.5 : max(kenya.geo$LATNUM) can't capture part of the North-West upper tip of the Kenya map
loc.pred = cbind(rep(xx, each = length(yy)), rep(yy, length(xx)))


#Removing the points (that are outside the Kenya) from the grid
#Import Kenya shape file
kenya_shape=readOGR('simulation_uniform/kenya_administrative1/ken_admbnda_adm1_iebc_20191031.shp')

#Convert prediction grid into a SpatialPointsDataFrame object 
grid=data.frame(xCoor = loc.pred[,1], yCoor = loc.pred[, 2])
grid <- SpatialPointsDataFrame(grid, data.frame(id=1:2500)) #we have 2500 points in the initial grid

#Remove the locations that are out of Kenya
grid=erase.point(grid, kenya_shape, inside = FALSE)

#Make a new loc.pred matrix from the remaining points
loc.pred = cbind(grid@coords[ ,1], grid@coords[ ,2])
nPred = dim(loc.pred)[1]


# Observation locations that are randomly chosen from a uniform distribution
loc.obs=data.frame(xCoor=runif(160, min=min(kenya.geo$LONGNUM), max=max(kenya.geo$LONGNUM)), yCoor=runif(160, min=min(kenya.geo$LATNUM), max=max(kenya.geo$LATNUM)))

#Removing the points (that are outside the Kenya) from the grid

#Convert observation locations into a SpatialPointsDataFrame object 
loc.obs=data.frame(xCoor = loc.obs[,1], yCoor = loc.obs[, 2])
loc.obs <- SpatialPointsDataFrame(loc.obs, data.frame(id=1:length(loc.obs[,1]))) #we have 2500 points in the initial grid

#Remove the locations that are out of Kenya
loc.obs=erase.point(loc.obs, kenya_shape, inside = FALSE)

#Make a new loc.pred matrix from the remaining points
loc.obs = cbind(loc.obs@coords[ ,1], loc.obs@coords[ ,2])
nLoc = dim(loc.obs)[1]

## Simulate spatial effect
# Desired parameters
range.sim = 1.5
sigma.sim = 1

# Covariance matrix
covMat = covFun(dMat = as.matrix(dist(rbind(loc.obs, loc.pred))),
                range = range.sim,
                stdDev = sigma.sim)

# Simulate spatial effect
L = t(chol(covMat))
u.sim = L%*%rnorm(dim(L)[1])

# Visual inspection (prediction grid)
df = data.frame(xCor = loc.pred[,1], yCor = loc.pred[, 2], u = u.sim[-(1:nLoc)])
df2 = data.frame(xCor = loc.obs[,1], yCor = loc.obs[, 2], u = u.sim[1:nLoc])
# ggplot() + geom_tile(data = data.frame(df, x1 = u.sim[-(1:nLoc)]), 
#                      aes(xCor, yCor, fill = x1)) +
#   coord_equal(xlim = c(33.5,42), ylim =c(-5,5))  +
#   geom_point(data = df2, aes(xCor,yCor), size = 0.2)
# 
# # Visual inspection (spatial observations)
# ggplot(data = df2) + geom_point(aes(xCor,yCor,color = u)) + coord_equal(xlim = c(33.5,42), ylim =c(-5,5))

# Save truth
myData = list(obs = df2,
              pred = df,
              truePar = list(range = range.sim,
                             spSD  = sigma.sim))

# Simulate observation model
# Desired parameters
sigma.nugget = sqrt(0.1)

# Data size
N = length(myData$obs$xCor)

# Covariate values
x = runif(N,-2,2)

# Coefficients of the fixed effects (intercept and slope)
alpha0 = 1
beta0 = 1

# Model Components
alpha   = rep(alpha0, N)
beta    = rep(beta0, N)
u       = myData$obs$u
epsilon = rnorm(N, 0, sigma.nugget)

# Assemble
lin.pred <- alpha + beta*x + u
y = lin.pred + epsilon
myData$obs$y = y
myData$obs$lin.pred = lin.pred
myData$obs$epsilon = epsilon
myData$truePar$epsSD = sigma.nugget
myData$obs$covar = matrix(x, nrow = length(x), ncol = 1)
myData$truePar$alpha = alpha
myData$truePar$beta  = beta

#Save  the simulated data set
save(myData, u.sim, kenya.geo, file="simulation_uniform/myDataOriginalUniform.RData")

# Visual inspection (observations)
ggplot(data = myData$obs) + geom_point(aes(xCor,yCor,color = y)) + coord_equal(xlim = c(33.5,42), ylim =c(-5,5))

## Analyse data with INLA
# Create mesh for inference
mesh.inla <- inla.mesh.2d(loc.domain = loc.pred,
                          max.edge = 0.13,
                          offset = -0.1)

# Create SPDE object
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.5)
spde.inla = inla.spde2.pcmatern(mesh = mesh.inla,
                                prior.range = prior.range,
                                prior.sigma = prior.sigma)

# Create observation matrix
A.obs = inla.spde.make.A(mesh = mesh.inla, loc = cbind(myData$obs$xCor, myData$obs$yCor))
A.pred = inla.spde.make.A(mesh = mesh.inla, loc = cbind(myData$pred$xCor, myData$pred$yCor))

# Create stack
stk.e <- inla.stack(tag='est',
                    data=list(y = as.vector(myData$obs$y)),
                    A=list(1, A.obs), 
                    effects=list(
                      list(intercept=rep(1, nLoc),
                           x=myData$obs$covar), 
                      list(idx.space = 1:spde.inla$n.spde)))
stk.pred = inla.stack(tag = 'pred',
                      data = list(y = rep(NA, nPred)),
                      A = list(1, A.pred),
                      effects = list(
                        list(intercept = rep(NA, nPred), 
                             x = rep(NA, nPred)),
                        list(idx.space = 1:spde.inla$n.spde)))
stk.full = inla.stack(stk.e, stk.pred)

# Formula
formula = y ~ intercept + x + f(idx.space, model = spde.inla) - 1

# Prior for nugget
prior.nugget = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.5), initial = log(1/0.5^2)))

# prior for fixed effects
prec.int = 1e-4
prec.beta = 1e-4

# Run INLA
# res.inla <- inla(formula = formula,
#                  data=inla.stack.data(stk.full, spde.inla = spde.inla),
#                  control.inla = list(int.strategy = 'grid',
#                                      dz = 0.75, 
#                                      diff.logdens = 15),
#                  control.family = list(hyper = prior.nugget),
#                  control.fixed = list(prec = prec.beta,
#                                       prec.intercept = prec.int),
#                  control.predictor=list(compute=TRUE,
#                                         A=inla.stack.A(stk.full)),
#                  verbose = TRUE)


res.inla <- inla(formula = formula,
                 data=inla.stack.data(stk.full, spde.inla = spde.inla),
                 control.inla = list(int.strategy = 'grid',
                                     dz = 0.75, 
                                     diff.logdens = 6),
                 control.family = list(hyper = prior.nugget),
                 control.fixed = list(prec = prec.beta,
                                      prec.intercept = prec.int),
                 control.predictor=list(compute=TRUE,
                                        A=inla.stack.A(stk.full)),
                 verbose = TRUE)



# Run inla hyperpar
res.inla.hyper = inla.hyperpar(res.inla, verbose = TRUE)

#Save results from INLA
save(res.inla, res.inla.hyper, stk.e, stk.pred, stk.full, file="simulation_uniform/INLA_originalUniform.RData")

## Analyse data with new STAN file
# Settings
nSamples = 80000

# Compile STAN code
model_stan = stan_model("simulation_uniform/stan_new.stan")

#creating the variable called "type". it refers to the location type (urban or rural)
type=rep(0, 100)
for (i in 1:100){
  if (kenya.geo@data[["URBAN_RURA"]][[i]]=="U") {
    type[[i]]=1 
  } else {
    type[[i]]=2
  }
}

# Create Stan data object
data.stan = list(N = length(myData$obs$y),
                     y = as.vector(myData$obs$y),
                     x = as.vector(myData$obs$covar),
                     L1 = myData$obs$xCor,
                     L2 = myData$obs$yCor,
                     priorSigma = prior.sigma,
                     priorRange = prior.range,
                     priorNugget = prior.nugget$prec$param,
                     nu = 1,
                     type = type,
                     max_urban = 6/111,
                     max_rural = 15/111,
                     pi = pi)

# Run Stan
res_stan = sampling(model_stan,
                        data = data.stan,
                        chains=1,
                        iter = nSamples,
                        init = 0.5,
                        thin = 8)

saveRDS(res_stan, file = "simulation_uniform/STAN_new_originalUniform.RDS")

## Sample spatial effects with new STAN results
# Extract parameter samples
thetaSample = extract(res_stan, pars = c("sdSpatial", "range", "sdNugget", "xCoor_new", "yCoor_new"))
sdSpatialSample = thetaSample[[1]]
rangeSample = thetaSample[[2]]
sdNuggetSample = thetaSample[[3]]
xCoor_new=thetaSample[[4]]
yCoor_new=thetaSample[[5]]
etaSample   = extract(res_stan, pars = c("eta"))[[1]]

# Initialize storage
uSample = matrix(NA, nrow = length(sdSpatialSample), ncol = nPred)

# Run through samples
for(i in 1:dim(uSample)[1]){
  print(i)
  # Extract current parameters
  range = rangeSample[i]
  sdSpatial = sdSpatialSample[i]
  sdNugget = sdNuggetSample[i]
 
  loc.obs = cbind(xCoor_new[i,], yCoor_new[i,])
  # Compute matrices
  dMat = as.matrix(dist(rbind(loc.obs, loc.pred)))
  SigAA = covFun(dMat[1:nLoc, 1:nLoc], range, sdSpatial) + sdNugget^2*diag(nrow = nLoc, ncol = nLoc)
  SigAB = covFun(dMat[1:nLoc, -(1:nLoc)], range, sdSpatial)
  SigBA = t(SigAB)
  SigBB = covFun(dMat[-(1:nLoc), -(1:nLoc)], range, sdSpatial)
  
  # Compute conditional distributions
  muC = matrix(0, nrow = nPred, ncol = 1) + SigBA%*%solve(SigAA, data.stan$y-etaSample[i,])
  SigC = SigBB - SigBA%*%solve(SigAA, SigAB)
  
  # Sample spatial effects at prediction locations
  L = t(chol(SigC))
  uSample[i,] = muC + as.vector(L%*%rnorm(nPred))
}

#Save results from analysis with new STAN file
save(res_stan, data.stan, uSample, thetaSample, sdSpatialSample, rangeSample, sdNuggetSample, etaSample, file="simulation_uniform/STAN_new_originalUniform.RData")

#Remove data from previous computation
rm(res_stan, data.stan)

## Analyse data with old STAN file
# Settings
nSamples = 80000

# Compile STAN code
model_stan = stan_model("simulation_uniform/stan_old.stan")

# Create Stan data object
data.stan = list(N = length(myData$obs$y),
                     y = as.vector(myData$obs$y),
                     x = as.vector(myData$obs$covar),
                     L1 = myData$obs$xCor,
                     L2 = myData$obs$yCor,
                     priorSigma = prior.sigma,
                     priorRange = prior.range,
                     priorNugget = prior.nugget$prec$param,
                     nu = 1)
# Run Stan
res_stan = sampling(model_stan,
                        data = data.stan,
                        chains=1,
                        iter = nSamples,
                        init = 0.5,
                        thin = 8)

saveRDS(res_stan, file = "simulation_uniform/STAN_old_originalUniform.RDS")

## Sample spatial effects with old STAN results
## Sample spatial effects
# Extract parameter samples
thetaSample = extract(res_stan, pars = c("sdSpatial", "range", "sdNugget"))
sdSpatialSample = thetaSample[[1]]
rangeSample = thetaSample[[2]]
sdNuggetSample = thetaSample[[3]]
etaSample   = extract(res_stan, pars = c("eta"))[[1]]

# Initialize storage
uSample = matrix(NA, nrow = length(sdSpatialSample), ncol = nPred)

# Run through samples
for(i in 1:dim(uSample)[1]){
  print(i)
  # Extract current parameters
  range = rangeSample[i]
  sdSpatial = sdSpatialSample[i]
  sdNugget = sdNuggetSample[i]
  
  # Compute matrices
  dMat = as.matrix(dist(rbind(loc.obs, loc.pred)))
  SigAA = covFun(dMat[1:nLoc, 1:nLoc], range, sdSpatial) + sdNugget^2*diag(nrow = nLoc, ncol = nLoc)
  SigAB = covFun(dMat[1:nLoc, -(1:nLoc)], range, sdSpatial)
  SigBA = t(SigAB)
  SigBB = covFun(dMat[-(1:nLoc), -(1:nLoc)], range, sdSpatial)
  
  # Compute conditional distributions
  muC = matrix(0, nrow = nPred, ncol = 1) + SigBA%*%solve(SigAA, data.stan$y-etaSample[i,])
  SigC = SigBB - SigBA%*%solve(SigAA, SigAB)
  
  # Sample spatial effects at prediction locations
  L = t(chol(SigC))
  uSample[i,] = muC + as.vector(L%*%rnorm(nPred))
}


#Save results from analysis with old STAN file
save(res_stan, data.stan, uSample, thetaSample, sdSpatialSample, rangeSample, sdNuggetSample, etaSample, file="simulation_uniform/STAN_old_originalUniform.RData")

quit()