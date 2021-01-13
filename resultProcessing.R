rm(list = ls())

# Load libraries
library(tidyverse)
library(rstan)
library(coda)
library(INLA)
library(rgdal)
library(ggplot2)
inla.setOption(pardiso.license = "simulation/pardiso.lic")


#POSTERIOR DISTRIBUTION GRAPHS WITH THE NEW STAN BASED ON THE ORIGINAL COORDINATES

#Set the working directory
directory="~/Desktop/STAN sampling"
setwd(directory)

#Results from new STAN with original coordinates
load("STAN_new_original.RData")

#Data set of Locations
load("myDataOriginal.RData")

#Data set of INLA Results
load("INLA_original.RData")

#Creating the graphs for comparing the posterior distributions of parameters
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.5)
prior.nugget = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.5), initial = log(1/0.5^2)))

#Trace plots of STAN runs
stan_trace(res_stan)

## Compare INLA and STAN New
# Intercept
marg.inla = res.inla.hyper$marginals.fixed[[1]]
samples.stan = extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())


#POSTERIOR DISTRIBUTION GRAPHS WITH THE OLD STAN BASED ON THE ORIGINAL COORDINATES

#Set the working directory
directory="~/Desktop/STAN sampling"
setwd(directory)

#Results from new STAN with original coordinates
load("STAN_old_original.RData")

#Data set of Locations
load("myDataOriginal.RData")

#Data set of INLA Results
load("INLA_original.RData")

#Creating the graphs for comparing the posterior distributions of parameters
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.5)
prior.nugget = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.5), initial = log(1/0.5^2)))

#Trace plots of STAN runs
stan_trace(res_stan)

## Compare INLA and STAN New
# Intercept
marg.inla = res.inla.hyper$marginals.fixed[[1]]
samples.stan = extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())


#POSTERIOR DISTRIBUTION GRAPHS WITH THE NEW STAN BASED ON THE JITTERED COORDINATES

#Set the working directory
directory="~/Desktop/STAN sampling"
setwd(directory)

#Results from new STAN with original coordinates
load("STAN_new_jittered.RData")

#Data set of Locations
load("myDataJittered.RData")

#Data set of INLA Results
load("INLA_jittered.RData")

#Creating the graphs for comparing the posterior distributions of parameters
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.5)
prior.nugget = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.5), initial = log(1/0.5^2)))

#Trace plots of STAN runs
stan_trace(res_stan)

## Compare INLA and STAN New
# Intercept
marg.inla = res.inla.hyper$marginals.fixed[[1]]
samples.stan = extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())


#POSTERIOR DISTRIBUTION GRAPHS WITH THE OLD STAN BASED ON THE JITTERED COORDINATES

#Set the working directory
directory="~/Desktop/STAN sampling"
setwd(directory)

#Results from new STAN with original coordinates
load("STAN_old_jittered.RData")

#Data set of Locations
load("myDataJittered.RData")

#Data set of INLA Results
load("INLA_jittered.RData")

#Creating the graphs for comparing the posterior distributions of parameters
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.5)
prior.nugget = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.5), initial = log(1/0.5^2)))

#Trace plots of STAN runs
stan_trace(res_stan)

## Compare INLA and STAN New
# Intercept
marg.inla = res.inla.hyper$marginals.fixed[[1]]
samples.stan = extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())


# 
# #PLOTTING THE SAMPLED LOCATIONS
# directory="~/Desktop/STAN sampling"
# setwd(directory)
# 
# load("STAN_new_original.RData")
# res_stan_original=res_stan
# rm(res_stan)
# 
# load("STAN_new_jittered.RData")
# res_stan_jittered=res_stan
# rm(res_stan)
# 
# #Data set of Locations
# load("myDataJittered.RData")
# 
# thetaSample_original = extract(res_stan_original)
# 
# sampledLoc_original1=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,1], Latitude = thetaSample_original[["yCoor_new"]][,1])
# sampledLoc_original2=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,2], Latitude = thetaSample_original[["yCoor_new"]][,2])
# sampledLoc_original3=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,3], Latitude = thetaSample_original[["yCoor_new"]][,3])
# sampledLoc_original4=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,4], Latitude = thetaSample_original[["yCoor_new"]][,4])
# sampledLoc_original5=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,5], Latitude = thetaSample_original[["yCoor_new"]][,5])
# sampledLoc_original6=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,6], Latitude = thetaSample_original[["yCoor_new"]][,6])
# sampledLoc_original7=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,7], Latitude = thetaSample_original[["yCoor_new"]][,7])
# sampledLoc_original8=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,8], Latitude = thetaSample_original[["yCoor_new"]][,8])
# sampledLoc_original9=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,9], Latitude = thetaSample_original[["yCoor_new"]][,9])
# sampledLoc_original10=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,10], Latitude = thetaSample_original[["yCoor_new"]][,10])
# 
# 
# thetaSample_jittered = extract(res_stan_jittered)
# 
# sampledLoc_jittered1=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,1], Latitude = thetaSample_jittered[["yCoor_new"]][,1])
# sampledLoc_jittered2=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,2], Latitude = thetaSample_jittered[["yCoor_new"]][,2])
# sampledLoc_jittered3=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,3], Latitude = thetaSample_jittered[["yCoor_new"]][,3])
# sampledLoc_jittered4=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,4], Latitude = thetaSample_jittered[["yCoor_new"]][,4])
# sampledLoc_jittered5=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,5], Latitude = thetaSample_jittered[["yCoor_new"]][,5])
# sampledLoc_jittered6=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,6], Latitude = thetaSample_jittered[["yCoor_new"]][,6])
# sampledLoc_jittered7=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,7], Latitude = thetaSample_jittered[["yCoor_new"]][,7])
# sampledLoc_jittered8=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,8], Latitude = thetaSample_jittered[["yCoor_new"]][,8])
# sampledLoc_jittered9=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,9], Latitude = thetaSample_jittered[["yCoor_new"]][,9])
# sampledLoc_jittered10=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,10], Latitude = thetaSample_jittered[["yCoor_new"]][,10])
# 
# #Extract the corresponding observed locations
# observed1=data.frame(Longitude = myData[["obs"]][["xCor"]][1], Latitude = myData[["obs"]][["yCor"]][1])
# observed2=data.frame(Longitude = myData[["obs"]][["xCor"]][2], Latitude = myData[["obs"]][["yCor"]][2])
# observed3=data.frame(Longitude = myData[["obs"]][["xCor"]][3], Latitude = myData[["obs"]][["yCor"]][3])
# observed4=data.frame(Longitude = myData[["obs"]][["xCor"]][4], Latitude = myData[["obs"]][["yCor"]][4])
# observed5=data.frame(Longitude = myData[["obs"]][["xCor"]][5], Latitude = myData[["obs"]][["yCor"]][5])
# observed6=data.frame(Longitude = myData[["obs"]][["xCor"]][6], Latitude = myData[["obs"]][["yCor"]][6])
# observed7=data.frame(Longitude = myData[["obs"]][["xCor"]][7], Latitude = myData[["obs"]][["yCor"]][7])
# observed8=data.frame(Longitude = myData[["obs"]][["xCor"]][8], Latitude = myData[["obs"]][["yCor"]][8])
# observed9=data.frame(Longitude = myData[["obs"]][["xCor"]][9], Latitude = myData[["obs"]][["yCor"]][9])
# observed10=data.frame(Longitude = myData[["obs"]][["xCor"]][10], Latitude = myData[["obs"]][["yCor"]][10])
# 
# 
# #Extract the corresponding jittered locations
# jittered1=data.frame(Longitude = myData[["jitt"]][["xCor"]][1], Latitude = myData[["jitt"]][["yCor"]][1])
# jittered2=data.frame(Longitude = myData[["jitt"]][["xCor"]][2], Latitude = myData[["jitt"]][["yCor"]][2])
# jittered3=data.frame(Longitude = myData[["jitt"]][["xCor"]][3], Latitude = myData[["jitt"]][["yCor"]][3])
# jittered4=data.frame(Longitude = myData[["jitt"]][["xCor"]][4], Latitude = myData[["jitt"]][["yCor"]][4])
# jittered5=data.frame(Longitude = myData[["jitt"]][["xCor"]][5], Latitude = myData[["jitt"]][["yCor"]][5])
# jittered6=data.frame(Longitude = myData[["jitt"]][["xCor"]][6], Latitude = myData[["jitt"]][["yCor"]][6])
# jittered7=data.frame(Longitude = myData[["jitt"]][["xCor"]][7], Latitude = myData[["jitt"]][["yCor"]][7])
# jittered8=data.frame(Longitude = myData[["jitt"]][["xCor"]][8], Latitude = myData[["jitt"]][["yCor"]][8])
# jittered9=data.frame(Longitude = myData[["jitt"]][["xCor"]][9], Latitude = myData[["jitt"]][["yCor"]][9])
# jittered10=data.frame(Longitude = myData[["jitt"]][["xCor"]][10], Latitude = myData[["jitt"]][["yCor"]][10])
# 
# 
# #Location1
# ggplot() + 
#   geom_point(data = sampledLoc_original1, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered1, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed1, aes(x = Longitude, y = Latitude, color = "Observed_Point"), size = 1) +
#   geom_point(data = jittered1, aes(x = Longitude, y = Latitude, color = "Jittered_Point"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","red", "yellow"))
# 
# #Location2
# ggplot() + 
#   geom_point(data = sampledLoc_original2, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered2, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed2, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location3
# ggplot() + 
#   geom_point(data = sampledLoc_original3, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered3, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed3, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location4
# ggplot() + 
#   geom_point(data = sampledLoc_original4, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered4, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed4, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location5
# ggplot() + 
#   geom_point(data = sampledLoc_original5, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered5, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed5, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location6
# ggplot() + 
#   geom_point(data = sampledLoc_original6, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered6, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed6, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location7
# ggplot() + 
#   geom_point(data = sampledLoc_original7, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered7, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed7, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location8
# ggplot() + 
#   geom_point(data = sampledLoc_original8, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered8, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed8, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location9
# ggplot() + 
#   geom_point(data = sampledLoc_original9, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered9, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed9, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
# 
# #Location10
# ggplot() + 
#   geom_point(data = sampledLoc_original10, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered10, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
#   geom_point(data = observed10, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
#   xlab('Longitude') +
#   ylab('Latitude')+
#   scale_colour_manual(values=c("orange","blue","yellow"))
#   
# 
# # contour lines
# 
# 
# p=ggplot() + 
#   geom_point(data = sampledLoc_original10, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
#   geom_point(data = sampledLoc_jittered10, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) 
# 
# 
# p+stat_density_2d_filled(alpha = 0.5)
# p + stat_density_2d(size = 0.25, colour = "black")
# 
# 
# 
#   xlab('Longitude') +
#   ylab('Latitude')
# geom_point(data = observed10, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1)
#   
# Var1=Longitude
# Var2=Latitude
# v3 = rbind(sampledLoc_original10, sampledLoc_jittered10)
# 
# library(ggplot2)
# 
# ggplot(v3, aes(x=Longitude, y=Latitude, colour="yellow")) +
#   stat_contour(binwidth=10) +
#   theme(panel.background=element_rect(fill="grey90")) +
#   theme(panel.grid=element_blank()) +
#   labs(title="Plot 1")
# 
# p2 = ggplot(v3, aes(x=Var1, y=Var2, z=value, colour=category)) +
#   stat_contour(aes(alpha=..level..), binwidth=10) +
#   theme(panel.background=element_rect(fill="white")) +
#   theme(panel.grid=element_blank()) +
#   labs(title="Plot 2")
# 
# p3 = ggplot(v3, aes(x=Var1, y=Var2, z=value, group=category)) +
#   stat_contour(aes(color=..level..), binwidth=10) +
#   scale_colour_gradient(low="white", high="#A1CD3A") +
#   theme(panel.background=element_rect(fill="grey50")) +
#   theme(panel.grid=element_blank()) +
#   labs(title="Plot 3")
# 
# p4 = ggplot(v3, aes(x=Var1, y=Var2, z=value, linetype=category)) +
#   stat_contour(aes(color=..level..), binwidth=10) +
#   scale_colour_gradient(low="white", high="#A1CD3A") +
#   theme(panel.background=element_rect(fill="grey50")) +
#   theme(panel.grid=element_blank()) +
#   labs(title="Plot 4")
# 
# library(gridExtra)
# ggsave(filename="plots.png", height=8, width=10,
#        plot=arrangeGrob(p1, p2, p3, p4, nrow=2, ncol=2))
# 

#Dawid Sebastiani Scores and RMSE for Old STAN Script based on Original Coordinates
rm(list = ls())
directory="~/Desktop/STAN sampling"
setwd(directory)

load("INLA_original.RData")
load("STAN_old_original.RData")
load("myDataOriginal.RData")

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaOrig=list()
DS_inlaOrig=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_original=(sum(unlist(sq_difference))/2500)^0.5

DS_stanOldOrig=list()
DS_stanOldOrig=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stanold_original=(sum(unlist(sq_difference))/2500)^0.5

boxplot(unlist(DS_stanOldOrig), horizontal=FALSE, ylim = c((min(c(unlist(DS_stanOldOrig),unlist(DS_inlaOrig)))-0.000000002), (max(c(unlist(DS_stanOldOrig), unlist(DS_inlaOrig)))+0.000000002)))
title(sub ="DS Scores STAN OLD", line = 0)
abline(h=mean(unlist(DS_stanOldOrig)), col ="red")   #its own mean      
abline(h=mean(unlist(DS_inlaOrig)), col ="blue") #mean crps obtained from INLA

save(rmse_stanold_original, file="results1.RData")

#Dawid Sebastiani Scores and RMSE for New STAN Script based on Original Coordinates
rm(list = ls())
directory="~/Desktop/STAN sampling"
setwd(directory)

load("INLA_original.RData")
load("STAN_new_original.RData")
load("myDataOriginal.RData")

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaOrig=list()
DS_inlaOrig=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_original=(sum(unlist(sq_difference))/2500)^0.5

DS_stanNewOrig=list()
DS_stanNewOrig=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stannew_original=(sum(unlist(sq_difference))/2500)^0.5
  
boxplot(unlist(DS_stanNewOrig),horizontal=FALSE, ylim = c((min(c(unlist(DS_stanNewOrig), DS_inlaOrig))-0.000000002), (max(c(unlist(DS_stanNewOrig), DS_inlaOrig))+0.000000002)))
title(sub ="DS Scores STAN NEW", line = 0)
abline(h=mean(DS_stanNewOrig), col ="red")   #its own mean      
abline(h=mean(DS_inlaOrig), col ="blue") #mean crps obtained from INLA
  
save(rmse_inla_original, rmse_stannew_original, file="results2.RData")

#Dawid Sebastiani Scores and RMSE for Old STAN Script based on Jittered Coordinates
rm(list = ls())
directory="~/Desktop/STAN sampling"
setwd(directory)

load("INLA_jittered.RData")
load("STAN_old_jittered.RData")
load("myDataJittered.RData")

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaJitt=list()
DS_inlaJitt=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_jittered=(sum(unlist(sq_difference))/2500)^0.5

DS_stanOldJitt=list()
DS_stanOldJitt=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stanold_jittered=(sum(unlist(sq_difference))/2500)^0.5

boxplot(unlist(DS_stanOldJitt),horizontal=FALSE, ylim = c((min(c(unlist(DS_stanOldJitt), DS_inlaJitt))-0.000000002), (max(c(unlist(DS_stanOldJitt), DS_inlaJitt))+0.000000002)))
title(sub ="DS Scores STAN OLD", line = 0)
abline(h=mean(DS_stanOldJitt), col ="red")   #its own mean      
abline(h=mean(DS_inlaJitt), col ="blue") #mean crps obtained from INLA

save(rmse_inla_jittered, rmse_stanold_jittered, file="results3.RData")

#Dawid Sebastiani Scores and RMSE for New STAN Script based on Jittered Coordinates
rm(list = ls())
directory="~/Desktop/STAN sampling"
setwd(directory)

load("INLA_jittered.RData")
load("STAN_new_jittered.RData")
load("myDataJittered.RData")

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaJitt=list()
DS_inlaJitt=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_jittered=(sum(unlist(sq_difference))/2500)^0.5

DS_stanNewJitt=list()
DS_stanNewJitt=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stannew_jittered=(sum(unlist(sq_difference))/2500)^0.5

boxplot(unlist(DS_stanNewJitt),horizontal=FALSE, ylim = c((min(c(unlist(DS_stanNewJitt), DS_inlaJitt))-0.000000002), (max(c(unlist(DS_stanNewJitt), DS_inlaJitt))+0.000000002)))
title(sub ="DS Scores STAN NEW", line = 0)
abline(h=mean(DS_stanNewJitt), col ="red")   #its own mean      
abline(h=mean(DS_inlaJitt), col ="blue") #mean crps obtained from INLA

save(rmse_stannew_jittered, file="results4.RData")

#Tabulation of RMSE values
rm(list = ls())
directory="~/Desktop/STAN sampling"
setwd(directory)

load("results1.RData")
load("results2.RData")
load("results3.RData")
load("results4.RData")

library(xtable)
rmse=data.frame(Type=c("Original", "Jittered"), INLA=c(rmse_inla_original, rmse_inla_jittered), STAN1=c(rmse_stanold_original, rmse_stanold_jittered), STAN2=c(rmse_stannew_original, rmse_stannew_jittered))
xtable(rmse)
