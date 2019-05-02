####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
BirdDat <- read.csv(text=getURL(gitstring))
####################################################################################################################################

x = BirdDat$RedtailedHawk
n = length(x)
N = 10^6

S=10000
#### DIC Pois-Gamma ####
alpha = 0.01
beta = 0.01

post.gamma.sample = rgamma(S, alpha + sum(x), n+beta)

poissonLogLik = function(x,theta){
  # total number of observations
  n = length(x)
  # equation
  logLik = (sum(x)*log(theta)-n*theta - sum(lfactorial(x) ))
  return(logLik)
}

hlp<-NULL
for(t in 1:S){
  hlp[t]<-poissonLogLik(x, post.gamma.sample[t])
}
lph<-poissonLogLik(x, mean(post.gamma.sample))
pdic=2*(lph-mean(hlp))
DIC_M1=-2*lph+2*pdic


#### DIC Beta ####
alpha = 1.019
beta = 1
post.beta.sample = rbeta(S, alpha + sum(x), n*N+ beta - sum(x))
binomLogLik = function(x, theta, N){
  n = length(x)
  logLik = sum(x) * log(theta) + (n*N - sum(x)) * log(1 - theta) + sum(lchoose(N, x))
  return(logLik)
}

hlp<-NULL
for(t in 1:S){
  hlp[t]<-binomLogLik(x, post.beta.sample[t],N)
}
lph<-binomLogLik(x, mean(post.beta.sample), N)
pdic=2*(lph-mean(hlp))
DIC_M2=-2*lph+2*pdic

(DIC = c(DIC_M1, DIC_M2))
abs(diff(DIC))


#### Gelfand and Ghosh Poisson ####
### Calculate Gelfand and Ghosh#####################
#first create predicted values for each of our posterior samples
#n columns (as many columns as we have data points) and S rows (number of posterior samples)
g.pred_values=matrix(0,S,n)
for(i in 1:length(post.gamma.sample)){
  #obtain one prediction for each hospital at each posterior sample
  g.pred_values[i,]=rpois(rep(1,n), rep(post.gamma.sample[i], length(x)))
}

#G term of gelfand and ghosh
g.g=sum((apply(g.pred_values, 2,mean)-x)^2)
g.p=sum(apply(g.pred_values,2,var))
gg_criterion_pois=g.g+g.p


#### Gelfand and Ghosh Binomial ####
### Calculate Gelfand and Ghosh#####################
#first create predicted values for each of our posterior samples
#n columns (as many columns as we have data points) and S rows (number of posterior samples)
b.pred_values=matrix(0,S,n)
for(i in 1:length(post.beta.sample)){
  #obtain one prediction for each hospital at each posterior sample
  b.pred_values[i,]=rbinom(rep(1,n), N, rep(post.beta.sample[i], length(x)))
}

#G term of gelfand and ghosh
b.g=sum((apply(b.pred_values, 2,mean)-x)^2)
b.p=sum(apply(b.pred_values,2,var))
gg_criterion_bin=b.g+b.p

(GGH = c(gg_criterion_pois, gg_criterion_bin))
abs(diff(GGH))


