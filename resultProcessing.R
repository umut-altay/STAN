rm(list = ls())

# Load libraries
library(rstan)
library(coda)
library(INLA)
library(rgdal)
library(ggplot2)
inla.setOption(pardiso.license = "simulation/pardiso.lic")

#Both R scripts (based on original and jittered locations) use same variable names
#But now the results need to be seperated to draw the graphs.
#import STAN results with original locations
setwd("~/Desktop")
load("STAN_new_original.RData")
#Save it with seperate name
res_stan_original=res_stan
save(res_stan_original, file="res_stan_original.RData")
rm(list = ls())


#import STAN results with jittered locations
load("STAN_new_jittered.RData")
#Save it with seperate name
res_stan_jittered=res_stan
save(res_stan_jittered, file="res_stan_jittered.RData")
rm(list = ls())

#Import the location (observation and prediction) data file and the both results
load("res_stan_original.RData")
load("res_stan_jittered.RData")
load("myDataOriginal.RData")

#Creating the graphs for comparing the posterior distributions of parameters
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.5)
prior.nugget = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.5), initial = log(1/0.5^2)))

#Trace plots of STAN runs
stan_trace(res_stan)


## Compare INLA and STAN
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

#Plotting the locations
#Extracting sampled coordinates based on the original locations
load("~/Desktop/res_stan_original.RData")

thetaSample_original = extract(res_stan_original)

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


#Extracting sampled coordinates based on the original locations
load("~/Desktop/res_stan_jittered.RData")
thetaSample_jittered = extract(res_stan_jittered)

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

#Extract the corresponding observed locations
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

#Location1
ggplot() + 
  geom_point(data = sampledLoc_original1, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered1, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed1, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location2
ggplot() + 
  geom_point(data = sampledLoc_original2, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered2, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed2, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location3
ggplot() + 
  geom_point(data = sampledLoc_original3, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered3, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed3, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location4
ggplot() + 
  geom_point(data = sampledLoc_original4, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered4, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed4, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location5
ggplot() + 
  geom_point(data = sampledLoc_original5, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered5, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed5, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location6
ggplot() + 
  geom_point(data = sampledLoc_original6, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered6, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed6, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location7
ggplot() + 
  geom_point(data = sampledLoc_original7, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered7, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed7, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location8
ggplot() + 
  geom_point(data = sampledLoc_original8, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered8, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed8, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location9
ggplot() + 
  geom_point(data = sampledLoc_original9, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered9, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed9, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

#Location10
ggplot() + 
  geom_point(data = sampledLoc_original10, aes(x = Longitude, y = Latitude, color = "Original"), size = 1) +
  geom_point(data = sampledLoc_jittered10, aes(x = Longitude, y = Latitude, color = "Jittered"), size = 1) +
  geom_point(data = observed10, aes(x = Longitude, y = Latitude, color = "Observed"), size = 1) +
  xlab('Longitude') +
  ylab('Latitude')+
  scale_colour_manual(values=c("orange","blue","yellow"))

grid=data.frame(Longitude = myData[["pred"]][["xCor"]], Latitude= myData[["pred"]][["yCor"]])

ggplot() + 
  geom_point(data = grid, aes(x = Longitude, y = Latitude), size = 1)


#Dawid Sebastiani Scores

load("~/Desktop/8 January Results/INLA_original.RData")
load("~/Desktop/8 January Results/STAN_old_original.RData")
load("~/Desktop/8 January Results/myDataOriginal.RData")

#Recreate INLA full stack
set.seed(2015104)

nLoc = 100
loc.pred = cbind(myData[["pred"]][["xCor"]], myData[["pred"]][["yCor"]])
nPred = dim(loc.pred)[1]

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


setwd("~/Desktop")
save(DS_inlaOrig, DS_stanOldOrig, rmse_inla_original, rmse_stanold_original, file="OriginalDSrmse_inla_stanOld.RData")
rm(list = ls())

load("~/Desktop/8 January Results/STAN_new_original.RData")
load("~/Desktop/8 January Results/myDataOriginal.RData")

DS_stanNewOrig=list()
DS_stanNewOrig=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stannew_original=(sum(unlist(sq_difference))/2500)^0.5

save(DS_stanNewOrig, rmse_stannew_original, file="OriginalDSrmse_stanNew.RData")
rm(list = ls())

#DS scores from jittered coordinates

load("~/Desktop/8 January Results/myDataJittered.RData")
load("~/Desktop/8 January Results/STAN_old_jittered.RData")
load("~/Desktop/8 January Results/INLA_jittered.RData")

#Recreate INLA full stack
set.seed(2015104)

nLoc = 100
loc.pred = cbind(myData[["pred"]][["xCor"]], myData[["pred"]][["yCor"]])
nPred = dim(loc.pred)[1]

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

#Extracting the Prediction Mean and Standard Deviation
index=inla.stack.index(stk.full, 'pred')$data
mean_pred = res.inla.hyper$summary.linear.predictor[index, "mean"]
sd_pred = res.inla.hyper$summary.linear.predictor[index, "sd"]

#DS scores and RMSE from original coordinates
DS_inlaJit=list()
DS_inlaJit=((myData[["pred"]][["u"]]-mean_pred)/sd_pred)^2+log(sd_pred^2)

sq_difference=list()
sq_difference=(mean_pred-myData[["pred"]][["u"]])^2
rmse_inla_jittered=(sum(unlist(sq_difference))/2500)^0.5


DS_stanOldJit=list()
DS_stanOldJit=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stanold_jittered=(sum(unlist(sq_difference))/2500)^0.5


setwd("~/Desktop")
save(DS_inlaJit, DS_stanOldJit, rmse_inla_jittered, rmse_stanold_jittered, file="JitteredDSrmse_inla_stanOld.RData")
rm(list = ls())

load("~/Desktop/8 January Results/STAN_new_jittered.RData")
load("~/Desktop/8 January Results/myDataJittered.RData")

DS_stanNewJit=list()
DS_stanNewJit=((myData[["pred"]][["u"]]-mean(uSample[1,]))/sd(uSample[1,]))^2+log(sd(uSample[1,])^2)

sq_difference=list()
sq_difference=(mean(uSample[1,])-myData[["pred"]][["u"]])^2
rmse_stannew_jittered=(sum(unlist(sq_difference))/2500)^0.5

save(DS_stanNewJit, rmse_stannew_jittered, file="JitteredDSrmse_stanNew.RData")
rm(list = ls())


#Graphs of DS Scores
#Original

load("~/Desktop/OriginalDSrmse_inla_stanOld.RData")
load("~/Desktop/OriginalDSrmse_stanNew.RData")

plot(DS_inlaOrig, DS_stanNewOrig, DS_stanOldOrig)


df_inla<- data.frame(DS=DS_inlaOrig,loc.index=1:2500)
df_stanOld<- data.frame(DS=DS_stanOldOrig,loc.index=1:2500)
df_stanNew<- data.frame(DS=DS_stanNewOrig,loc.index=1:2500)

ggplot() + 
  geom_point(data = df_inla, aes(x = loc.index, y = DS, color = "inla"), size = 1) +
  geom_point(data = df_stanOld, aes(x = loc.index, y = DS, color = "stanOld"), size = 1) +
  geom_point(data = df_stanNew, aes(x = loc.index, y = DS, color = "stanNew"), size = 1) +
  xlab('Prediction Location Number') +
  ylab('DS Score')+
  scale_colour_manual(values=c("red","yellow","blue"))
rm(list = ls())

#Jittered

load("~/Desktop/JitteredDSrmse_inla_stanOld.RData")
load("~/Desktop/JitteredDSrmse_stanNew.RData")

df_inla<- data.frame(DS=DS_inlaJit,loc.index=1:2500)
df_stanOld<- data.frame(DS=DS_stanOldJit,loc.index=1:2500)
df_stanNew<- data.frame(DS=DS_stanNewJit,loc.index=1:2500)

ggplot() + 
  geom_point(data = df_inla, aes(x = loc.index, y = DS, color = "inla"), size = 1) +
  geom_point(data = df_stanOld, aes(x = loc.index, y = DS, color = "stanOld"), size = 1) +
  geom_point(data = df_stanNew, aes(x = loc.index, y = DS, color = "stanNew"), size = 1) +
  xlab('Prediction Location Number') +
  ylab('DS Score')+
  scale_colour_manual(values=c("red","yellow","blue"))



#RMSE
load("~/Desktop/OriginalDSrmse_inla_stanOld.RData")
load("~/Desktop/OriginalDSrmse_stanNew.RData")
library(xtable)
rmse=data.frame(type=c("Original", "Jittered"), INLA=c(rmse_inla_original, rmse_inla_jittered), STAN_Old=c(rmse_stanold_original, rmse_stanold_jittered), STAN_New=c(rmse_stannew_original, rmse_stannew_jittered))
xtable(rmse)
