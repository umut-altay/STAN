rm(list = ls())

# Load libraries
#library(tidyverse)
library(rstan)
library(coda)
library(INLA)
library(rgdal)
library(ggplot2)
library(plyr)
inla.setOption(pardiso.license = "simulation/pardiso.lic")

#POSTERIOR DISTRIBUTION GRAPHS WITH THE NEW STAN BASED ON THE ORIGINAL COORDINATES

#Set the working directory
directory="~/Desktop/STAN sampling/UNIFORM"
#directory="~/Desktop/STAN sampling/GRID"
#directory="~/Desktop/STAN sampling"

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
samples.stan = rstan::extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = rstan::extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = rstan::extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = rstan::extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = rstan::extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())

#POSTERIOR DISTRIBUTION GRAPHS WITH THE OLD STAN BASED ON THE ORIGINAL COORDINATES

#Set the working directory
# directory="~/Desktop/STAN sampling"
# setwd(directory)

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
samples.stan = rstan::extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = rstan::extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = rstan::extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = rstan::extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = rstan::extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())

#POSTERIOR DISTRIBUTION GRAPHS WITH THE NEW STAN BASED ON THE JITTERED COORDINATES

#Set the working directory
# directory="~/Desktop/STAN sampling"
# setwd(directory)

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
samples.stan = rstan::extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = rstan::extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = rstan::extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = rstan::extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = rstan::extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())

#POSTERIOR DISTRIBUTION GRAPHS WITH THE OLD STAN BASED ON THE JITTERED COORDINATES

#Set the working directory
# directory="~/Desktop/STAN sampling"
# setwd(directory)

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
samples.stan = rstan::extract(res_stan, pars = "alpha")[[1]]
hist(samples.stan, 50, freq = F, main = "Intercept")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# Covariate
marg.inla = res.inla.hyper$marginals.fixed[[2]]
samples.stan = rstan::extract(res_stan, pars = "beta")[[1]]
hist(samples.stan, 50, freq = F, main = "Beta")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-1000, 1000, length.out = 10000)
yy = dnorm(xx, mean = 0, sd = 100)
lines(xx, yy, lwd = 2, col = "blue")

# log(Nugget std.dev.)
marg.inla = inla.tmarginal(function(x) {-x/2}, res.inla.hyper$internal.marginals.hyperpar[[1]])
samples.stan = rstan::extract(res_stan, pars = "theta[1]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Nugget)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.nugget$prec$param[2])/prior.nugget$prec$param[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Range)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[2]]
samples.stan = rstan::extract(res_stan, pars = "theta[3]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(Range)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.range[2])*prior.range[1]
yy = lam*exp(-lam*exp(-xx))*exp(-xx)
lines(xx, yy, lwd = 2, col = "blue")

# log(Spatial std.dev.)
marg.inla = res.inla.hyper$internal.marginals.hyperpar[[3]]
samples.stan = rstan::extract(res_stan, pars = "theta[2]")[[1]]
hist(samples.stan, 50, freq = F, main = "log(std.dev. Spatial)")
lines(marg.inla, lwd = 2, col = "red")
xx = seq(-20, 20, length.out = 1000)
lam = -log(prior.sigma[2])/prior.sigma[1]
yy = lam*exp(-lam*exp(xx))*exp(xx)
lines(xx, yy, lwd = 2, col = "blue")

rm(list = ls())

#PLOTTING THE SAMPLED LOCATIONS
rm(list = ls())

library(rstan)

# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("STAN_new_original.RData")
res_stan_original=res_stan
rm(res_stan)

load("STAN_new_jittered.RData")
res_stan_jittered=res_stan
rm(res_stan)

#Data set of Locations
load("myDataOriginal.RData")

#Extract data from STAN results
thetaSample_original = rstan::extract(res_stan_original)

sampledLoc_original1=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,1], Latitude = thetaSample_original[["yCoor_new"]][,1])
sampledLoc_original2=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,2], Latitude = thetaSample_original[["yCoor_new"]][,2])
sampledLoc_original3=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,3], Latitude = thetaSample_original[["yCoor_new"]][,3])
sampledLoc_original4=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,4], Latitude = thetaSample_original[["yCoor_new"]][,4])
sampledLoc_original5=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,5], Latitude = thetaSample_original[["yCoor_new"]][,5])
sampledLoc_original6=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,6], Latitude = thetaSample_original[["yCoor_new"]][,6])
sampledLoc_original7=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,7], Latitude = thetaSample_original[["yCoor_new"]][,7])
sampledLoc_original8=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,8], Latitude = thetaSample_original[["yCoor_new"]][,8])
sampledLoc_original9=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,9], Latitude = thetaSample_original[["yCoor_new"]][,9])
sampledLoc_original10=data.frame(Longitude = thetaSample_original[["xCoor_new"]][,10], Latitude = thetaSample_original[["yCoor_new"]][,10])

thetaSample_jittered = rstan::extract(res_stan_jittered)

sampledLoc_jittered1=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,1], Latitude = thetaSample_jittered[["yCoor_new"]][,1])
sampledLoc_jittered2=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,2], Latitude = thetaSample_jittered[["yCoor_new"]][,2])
sampledLoc_jittered3=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,3], Latitude = thetaSample_jittered[["yCoor_new"]][,3])
sampledLoc_jittered4=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,4], Latitude = thetaSample_jittered[["yCoor_new"]][,4])
sampledLoc_jittered5=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,5], Latitude = thetaSample_jittered[["yCoor_new"]][,5])
sampledLoc_jittered6=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,6], Latitude = thetaSample_jittered[["yCoor_new"]][,6])
sampledLoc_jittered7=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,7], Latitude = thetaSample_jittered[["yCoor_new"]][,7])
sampledLoc_jittered8=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,8], Latitude = thetaSample_jittered[["yCoor_new"]][,8])
sampledLoc_jittered9=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,9], Latitude = thetaSample_jittered[["yCoor_new"]][,9])
sampledLoc_jittered10=data.frame(Longitude = thetaSample_jittered[["xCoor_new"]][,10], Latitude = thetaSample_jittered[["yCoor_new"]][,10])

# #Extract the corresponding observed locations
observed1=data.frame(Longitude = myData[["obs"]][["xCor"]][1], Latitude = myData[["obs"]][["yCor"]][1])
observed2=data.frame(Longitude = myData[["obs"]][["xCor"]][2], Latitude = myData[["obs"]][["yCor"]][2])
observed3=data.frame(Longitude = myData[["obs"]][["xCor"]][3], Latitude = myData[["obs"]][["yCor"]][3])
observed4=data.frame(Longitude = myData[["obs"]][["xCor"]][4], Latitude = myData[["obs"]][["yCor"]][4])
observed5=data.frame(Longitude = myData[["obs"]][["xCor"]][5], Latitude = myData[["obs"]][["yCor"]][5])
observed6=data.frame(Longitude = myData[["obs"]][["xCor"]][6], Latitude = myData[["obs"]][["yCor"]][6])
observed7=data.frame(Longitude = myData[["obs"]][["xCor"]][7], Latitude = myData[["obs"]][["yCor"]][7])
observed8=data.frame(Longitude = myData[["obs"]][["xCor"]][8], Latitude = myData[["obs"]][["yCor"]][8])
observed9=data.frame(Longitude = myData[["obs"]][["xCor"]][9], Latitude = myData[["obs"]][["yCor"]][9])
observed10=data.frame(Longitude = myData[["obs"]][["xCor"]][10], Latitude = myData[["obs"]][["yCor"]][10])

#GGplot Figures
#NOTE: the longitude and latitude values in the geom_point(aes()) part at the end
#are directly taken from observed1, 2, 3,... data frame and copy pasted there.
#Because the ggplot code didn't recognize them when observed1[1,1] and observed1[1,2] is used instead

#Location1
X=c(sampledLoc_original1[,1], sampledLoc_jittered1[,1])
Y=c(sampledLoc_original1[,2], sampledLoc_jittered1[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][1], myData[["obs"]][["yCor"]][1]), color = "blue", size = 5) +
  coord_fixed(ratio = 1)


#Location2
X=c(sampledLoc_original2[,1], sampledLoc_jittered2[,1])
Y=c(sampledLoc_original2[,2], sampledLoc_jittered2[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][2], myData[["obs"]][["yCor"]][2]), color = "blue", size = 5) + 
  coord_fixed(ratio = 1)

#Location3
X=c(sampledLoc_original3[,1], sampledLoc_jittered3[,1])
Y=c(sampledLoc_original3[,2], sampledLoc_jittered3[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][3], myData[["obs"]][["yCor"]][3]), color = "blue", size = 5) +
  coord_fixed(ratio = 1) 

#Location4
X=c(sampledLoc_original4[,1], sampledLoc_jittered4[,1])
Y=c(sampledLoc_original4[,2], sampledLoc_jittered4[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][4], myData[["obs"]][["yCor"]][4]), color = "blue", size = 5) +
  coord_fixed(ratio = 1) 


#Location5
X=c(sampledLoc_original5[,1], sampledLoc_jittered5[,1])
Y=c(sampledLoc_original5[,2], sampledLoc_jittered5[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][5], myData[["obs"]][["yCor"]][5]), color = "blue", size = 5) +
  coord_fixed(ratio = 1)   

#Location6
X=c(sampledLoc_original6[,1], sampledLoc_jittered6[,1])
Y=c(sampledLoc_original6[,2], sampledLoc_jittered6[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][6], myData[["obs"]][["yCor"]][6]), color = "blue", size = 5) +
  coord_fixed(ratio = 1)   

#Location7
X=c(sampledLoc_original7[,1], sampledLoc_jittered7[,1])
Y=c(sampledLoc_original7[,2], sampledLoc_jittered7[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][7], myData[["obs"]][["yCor"]][7]), color = "blue", size = 5) +
  coord_fixed(ratio = 1)   

#Location8
X=c(sampledLoc_original8[,1], sampledLoc_jittered8[,1])
Y=c(sampledLoc_original8[,2], sampledLoc_jittered8[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][8], myData[["obs"]][["yCor"]][8]), color = "blue", size = 5) +
  coord_fixed(ratio = 1) 

#Location9
X=c(sampledLoc_original9[,1], sampledLoc_jittered9[,1])
Y=c(sampledLoc_original9[,2], sampledLoc_jittered9[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][9], myData[["obs"]][["yCor"]][9]), color = "blue", size = 5) +
  coord_fixed(ratio = 1) 

#Location10
X=c(sampledLoc_original10[,1], sampledLoc_jittered10[,1])
Y=c(sampledLoc_original10[,2], sampledLoc_jittered10[,2])
Label = c(rep('Original', 5000), rep('Jittered', 5000))
plot_data <-data.frame(X, Y, Label)
ggplot(plot_data, aes(X, Y, group = Label)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Label),
                  bins = 20) +
  geom_point(aes(myData[["obs"]][["xCor"]][10], myData[["obs"]][["yCor"]][10]), color = "blue", size = 5) +
  coord_fixed(ratio = 1) 
# #
# #


#Dawid Sebastiani Scores and RMSE for Old STAN Script based on Original Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("INLA_original.RData")
load("STAN_old_original.RData")
load("myDataOriginal.RData")

library(rstan)
library(INLA)
library(plyr) #for colwise(sd)(uSample) calculation
library(RANN)

#extract observation locations and prediction locations from myData
loc.obs = cbind(myData[["obs"]][["xCor"]], myData[["obs"]][["yCor"]])
loc.pred=cbind(myData[["pred"]][["xCor"]], myData[["pred"]][["yCor"]])
#Find nearest observation location for each prediction location (and the distances between them (as degrees))
nearestOriginal=as.data.frame(RANN::nn2(loc.obs[,c(1,2)],loc.pred[,c(1,2)],k=1))

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaOrig=list()
DS_inlaOrig=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_original=(sum(unlist(sq_difference))/1433)^0.5

sDev=apply(uSample, 2, sd)
DS_stanOldOrig=list()
DS_stanOldOrig=((myData[["pred"]][["u"]]-colMeans(uSample))/sDev)^2+log(sDev^2)

sq_difference=list()
sq_difference=(colMeans(uSample)-myData[["pred"]][["u"]])^2 
rmse_stanold_original=(sum(unlist(sq_difference))/1433)^0.5

#DS without separation with respect to neighbourhood distances
DS_inlaOrigAVG=mean(unlist(DS_inlaOrig))
DS_stanOldOrigAvg=mean(unlist(DS_stanOldOrig))

#Nearest neighbour results and DS scores becomes a single data frame
predPoint.index=c(1:length(loc.pred[,1]))
distClass_oldOriginal=data.frame(pred.point = predPoint.index, nnDistance = nearestOriginal[,2], DS_inla = DS_inlaOrig, DS_stan = DS_stanOldOrig)

boxplot(unlist(DS_stanOldOrig), horizontal=FALSE, ylim = c((min(c(unlist(DS_stanOldOrig),unlist(DS_inlaOrig)))-0.000000002), (max(c(unlist(DS_stanOldOrig), unlist(DS_inlaOrig)))+0.000000002)))
title(sub ="DS Scores STAN old-original", line = 0)
abline(h=mean(unlist(DS_stanOldOrig)), col ="red")   #its own mean      
abline(h=mean(unlist(DS_inlaOrig)), col ="blue") #mean crps obtained from INLA

save(rmse_stanold_original, distClass_oldOriginal, file="resultsOldOriginal.RData")
save(DS_inlaOrigAVG, DS_stanOldOrigAvg, DS_inlaOrig, DS_stanOldOrig, file = "DSoldOriginalUniform.RData")

#Dawid Sebastiani Scores and RMSE for New STAN Script based on Original Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("INLA_original.RData")
load("STAN_new_original.RData")
load("myDataOriginal.RData")

library(tidyverse)
library(rstan)
library(INLA)
library(plyr)
library(RANN)

#extract observation locations and prediction locations from myData
loc.obs = cbind(myData[["obs"]][["xCor"]], myData[["obs"]][["yCor"]])
loc.pred=cbind(myData[["pred"]][["xCor"]], myData[["pred"]][["yCor"]])
#Find nearest observation location for each prediction location (and the distances between them (as degrees))
nearestOriginal=as.data.frame(RANN::nn2(loc.obs[,c(1,2)],loc.pred[,c(1,2)],k=1))

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaOrig=list()
DS_inlaOrig=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_original=(sum(unlist(sq_difference))/1433)^0.5

sDev=apply(uSample, 2, sd)
DS_stanNewOrig=list()
DS_stanNewOrig=((myData[["pred"]][["u"]]-colMeans(uSample))/sDev)^2+log(sDev^2)


#DS without separation with respect to neighbourhood distances
DS_stanNewOrigAvg=mean(unlist(DS_stanNewOrig))


sq_difference=list()
sq_difference=(colMeans(uSample)-myData[["pred"]][["u"]])^2 
rmse_stannew_original=(sum(unlist(sq_difference))/1433)^0.5
  
#Nearest neighbour results and DS scores becomes a single data frame
predPoint.index=c(1:length(loc.pred[,1]))
distClass_newOriginal=data.frame(pred.point = predPoint.index, nnDistance = nearestOriginal[,2], DS_inla = DS_inlaOrig, DS_stan = DS_stanNewOrig)

boxplot(unlist(DS_stanNewOrig),horizontal=FALSE, ylim = c((min(c(unlist(DS_stanNewOrig), DS_inlaOrig))-0.000000002), (max(c(unlist(DS_stanNewOrig), DS_inlaOrig))+0.000000002)))
title(sub ="DS Scores STAN new-original", line = 0)
abline(h=mean(DS_stanNewOrig), col ="red")   #its own mean      
abline(h=mean(DS_inlaOrig), col ="blue") #mean crps obtained from INLA
  
save(DS_stanNewOrigAvg, rmse_inla_original, rmse_stannew_original, distClass_newOriginal, file="resultsNewOriginal.RData")
save(DS_stanNewOrig, file = "DSnewOriginalUniform.RData")

#Dawid Sebastiani Scores and RMSE for Old STAN Script based on Jittered Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("INLA_jittered.RData")
load("STAN_old_jittered.RData")
load("myDataOriginal.RData")

library(tidyverse)
library(rstan)
library(INLA)
library(plyr)
library(RANN)

#extract observation locations and prediction locations from myData
loc.obs = cbind(myData[["obs"]][["xCor"]], myData[["obs"]][["yCor"]])
loc.pred=cbind(myData[["pred"]][["xCor"]], myData[["pred"]][["yCor"]])
#Find nearest observation location for each prediction location (and the distances between them (as degrees))
nearestJittered=as.data.frame(RANN::nn2(loc.obs[,c(1,2)],loc.pred[,c(1,2)],k=1))

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaJitt=list()
DS_inlaJitt=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_jittered=(sum(unlist(sq_difference))/1433)^0.5

sDev=apply(uSample, 2, sd)
DS_stanOldJitt=list()
DS_stanOldJitt=((myData[["pred"]][["u"]]-colMeans(uSample))/sDev)^2+log(sDev^2)


#DS without separation with respect to neighbourhood distances
DS_inlaJittAVG=mean(unlist(DS_inlaJitt))
DS_stanOldJittAvg=mean(unlist(DS_stanOldJitt))



sq_difference=list()
sq_difference=(colMeans(uSample)-myData[["pred"]][["u"]])^2 
rmse_stanold_jittered=(sum(unlist(sq_difference))/1433)^0.5

#Nearest neighbour results and DS scores becomes a single data frame
predPoint.index=c(1:length(loc.pred[,1]))
distClass_oldJittered=data.frame(pred.point = predPoint.index, nnDistance = nearestJittered[,2], DS_inla = DS_inlaJitt, DS_stan = DS_stanOldJitt)

boxplot(unlist(DS_stanOldJitt),horizontal=FALSE, ylim = c((min(c(unlist(DS_stanOldJitt), DS_inlaJitt))-0.000000002), (max(c(unlist(DS_stanOldJitt), DS_inlaJitt))+0.000000002)))
title(sub ="DS Scores STAN old-jittered", line = 0)
abline(h=mean(DS_stanOldJitt), col ="red")   #its own mean      
abline(h=mean(DS_inlaJitt), col ="blue") #mean crps obtained from INLA

save(rmse_inla_jittered, rmse_stanold_jittered, distClass_oldJittered, file="resultsOldJittered.RData")
save(DS_inlaJittAVG, DS_stanOldJittAvg, DS_inlaJitt, DS_stanOldJitt, file = "DSoldJitteredUniform.RData")

#Dawid Sebastiani Scores and RMSE for New STAN Script based on Jittered Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("INLA_jittered.RData")
load("STAN_new_jittered.RData")
load("myDataOriginal.RData")

library(tidyverse)
library(rstan)
library(INLA)
library(plyr)
library(RANN)

#extract observation locations and prediction locations from myData
loc.obs = cbind(myData[["obs"]][["xCor"]], myData[["obs"]][["yCor"]])
loc.pred=cbind(myData[["pred"]][["xCor"]], myData[["pred"]][["yCor"]])
#Find nearest observation location for each prediction location (and the distances between them (as degrees))
nearestJittered=as.data.frame(RANN::nn2(loc.obs[,c(1,2)],loc.pred[,c(1,2)],k=1))

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaJitt=list()
DS_inlaJitt=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_jittered=(sum(unlist(sq_difference))/1433)^0.5

sDev=apply(uSample, 2, sd)
DS_stanNewJitt=list()
DS_stanNewJitt=((myData[["pred"]][["u"]]-colMeans(uSample))/sDev)^2+log(sDev^2)


#DS without separation with respect to neighbourhood distances
DS_stanNewJittAvg=mean(unlist(DS_stanNewJitt))


sq_difference=list()
sq_difference=(colMeans(uSample)-myData[["pred"]][["u"]])^2 
rmse_stannew_jittered=(sum(unlist(sq_difference))/1433)^0.5

#Nearest neighbour results and DS scores becomes a single data frame
predPoint.index=c(1:length(loc.pred[,1]))
distClass_newJittered=data.frame(pred.point = predPoint.index, nnDistance = nearestJittered[,2], DS_inla = DS_inlaJitt, DS_stan = DS_stanNewJitt)

boxplot(unlist(DS_stanNewJitt),horizontal=FALSE, ylim = c((min(c(unlist(DS_stanNewJitt), DS_inlaJitt))-0.000000002), (max(c(unlist(DS_stanNewJitt), DS_inlaJitt))+0.000000002)))
title(sub ="DS Scores STAN new-jittered", line = 0)
abline(h=mean(DS_stanNewJitt), col ="red")   #its own mean      
abline(h=mean(DS_inlaJitt), col ="blue") #mean crps obtained from INLA

save(rmse_stannew_jittered, distClass_newJittered, file="resultsNewJittered.RData")
save(DS_stanNewJitt, file = "DSnewJitteredUniform.RData")

#Tabulation of RMSE values
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("resultsOldOriginal.RData")
load("resultsNewOriginal.RData")
load("resultsOldJittered.RData")
load("resultsNewJittered.RData")

library(xtable)
rmse=data.frame(Type=c("Original", "Jittered"), INLA=c(rmse_inla_original, rmse_inla_jittered), STAN1=c(rmse_stanold_original, rmse_stanold_jittered), STAN2=c(rmse_stannew_original, rmse_stannew_jittered))
xtable(rmse)


#Splitting the data into 10 classes with respect to the distances and calculating average DS scores for each class

#For New STAN and Original Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("resultsOldOriginal.RData")
load("resultsNewOriginal.RData")
load("resultsOldJittered.RData")
load("resultsNewJittered.RData")

#Define the class size
#find a class size for 10 classes
size=(max(distClass_newOriginal$nnDistance)-min(distClass_newOriginal$nnDistance))/10
#breaks
b=c(-Inf, 1:9, Inf)
b=size*b
km=b*111
#labels
l=c(1:10) # labels should be 1 less than the breaks

#extract the nearest neighbour distance column form the data frame
nnd=distClass_newOriginal$nnDistance
cls=cut(nnd, breaks = b, labels = l)
distClass_newOriginal$cls=cls
avrScores_newOriginal=aggregate(distClass_newOriginal[, 3:4], list(distClass_newOriginal$cls), mean)
avrScores_newOriginal=data.frame(avrScores_newOriginal,Lower = km[1:10], Upper = km[2:11])
library(xtable)
xtable(avrScores_newOriginal)


#For Old STAN and Original Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("resultsOldOriginal.RData")
load("resultsNewOriginal.RData")
load("resultsOldJittered.RData")
load("resultsNewJittered.RData")

#Define the class size
#find a class size for 10 classes
size=(max(distClass_oldOriginal$nnDistance)-min(distClass_oldOriginal$nnDistance))/10
#breaks
b=c(-Inf, 1:9, Inf)
b=size*b
km=b*111
#labels
l=c(1:10) # labels should be 1 less than the breaks

#extract the nearest neighbour distance column form the data frame
nnd=distClass_oldOriginal$nnDistance
cls=cut(nnd, breaks = b, labels = l)
distClass_oldOriginal$cls=cls
avrScores_oldOriginal=aggregate(distClass_oldOriginal[, 3:4], list(distClass_oldOriginal$cls), mean)
avrScores_oldOriginal=data.frame(avrScores_oldOriginal,Lower = km[1:10], Upper = km[2:11])

library(xtable)
xtable(avrScores_oldOriginal)


#For New STAN and Jittered Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("resultsOldOriginal.RData")
load("resultsNewOriginal.RData")
load("resultsOldJittered.RData")
load("resultsNewJittered.RData")

#Define the class size
#find a class size for 10 classes
size=(max(distClass_newJittered$nnDistance)-min(distClass_newJittered$nnDistance))/10
#breaks
b=c(-Inf, 1:9, Inf)
b=size*b
km=b*111
#labels
l=c(1:10) # labels should be 1 less than the breaks

#extract the nearest neighbour distance column form the data frame
nnd=distClass_newJittered$nnDistance
cls=cut(nnd, breaks = b, labels = l)
distClass_newJittered$cls=cls
avrScores_newJittered=aggregate(distClass_newJittered[, 3:4], list(distClass_newJittered$cls), mean)
avrScores_newJittered=data.frame(avrScores_newJittered,Lower = km[1:10], Upper = km[2:11])

library(xtable)
xtable(avrScores_newJittered)


#For Old STAN and Jittered Coordinates
rm(list = ls())
# directory="~/Desktop/STAN sampling"
# setwd(directory)

load("resultsOldOriginal.RData")
load("resultsNewOriginal.RData")
load("resultsOldJittered.RData")
load("resultsNewJittered.RData")

#Define the class size
#find a class size for 10 classes
size=(max(distClass_oldJittered$nnDistance)-min(distClass_oldJittered$nnDistance))/10
#breaks
b=c(-Inf, 1:9, Inf)
b=size*b
km=b*111
#labels
l=c(1:10) # labels should be 1 less than the breaks

#extract the nearest neighbour distance column form the data frame
nnd=distClass_oldJittered$nnDistance
cls=cut(nnd, breaks = b, labels = l)
distClass_oldJittered$cls=cls
avrScores_oldJittered=aggregate(distClass_oldJittered[, 3:4], list(distClass_oldJittered$cls), mean)
avrScores_oldJittered=data.frame(avrScores_oldJittered,Lower = km[1:10], Upper = km[2:11])

library(xtable)
xtable(avrScores_oldJittered)



#Plotting Original and Jittered DS Scores Against Each Other

#Initial
directory="~/Desktop/STAN sampling"
setwd(directory)

load("DSnewJittered.RData")
load("DSoldJittered.RData")
load("DSnewOriginal.RData")
load("DSnewOriginal.RData")

#Uniform
directory="~/Desktop/STAN sampling/UNIFORM"
setwd(directory)

load("DSnewJitteredUniform.RData")
load("DSoldJitteredUniform.RData")
load("DSnewOriginalUniform.RData")
load("DSnewOriginalUniform.RData")

#Grid
directory="~/Desktop/STAN sampling/GRID"
setwd(directory)

load("DSnewJitteredGrid.RData")
load("DSoldJitteredGrid.RData")
load("DSnewOriginalGrid.RData")
load("DSnewOriginalGrid.RData")

plot(DS_stanOldOrig, DS_stanOldJitt)












