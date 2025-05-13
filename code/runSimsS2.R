set.seed(8)
options(warn=2)

source("Simulation_function.R")

#####################################################
# SET THE VALUES OF THE PARAMETERS THAT DONT CHANGE #
#####################################################


nSim<-250
nRand<-200  

nMeas<-4
nTrial<-8


beta.real<-c(0.4,-0.15)
sigmau<-1.5
sigmav<-0.3

z <- as.factor(1:nRand)
x0<-runif(nRand,0,1.5)
dat=data.frame(z,x0)
dat$x1<-dat$x0+1
dat$x2<-dat$x1+1
dat$x3<-dat$x2+3


unorm <- rnorm(nRand,0,sigmau) # Simulate the random effects
vnorm<-rnorm(nRand,0,sigmav) 

weibull<-c(0.1,1.6)
cens<-0.1

#############
# SIMULATE  #
#############

out <- simulation.function(nSim=nSim,nRand= nRand,nMeas=nMeas,nTrial= nTrial,beta= beta.real, unorm=unorm,vnorm=vnorm,phi=phi, sigmau=sigmau,sigmav=sigmav,weibull =  weibull,alpha=alpha,cens=cens,dat=dat) # get the results of the simulations for the given values


