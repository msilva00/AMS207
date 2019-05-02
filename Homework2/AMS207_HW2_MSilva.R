horsekick = c(3,	5,	7,	9,	10,
              18,	6,	14,	11,	9,
              5,	11,	15,	6,	11,
              17,	12,	15,	8,	4)
mean(horsekick)
alpha = 9.8
n = 20
beta = 1
x.grid = seq(6,14, length = 1000)
n=20
alpha = 9.8
beta = 1
plot(x.grid, dgamma(x.grid, sum(horsekick) + alpha, n + beta), type = 'l')
N = 10^6
colors = c("red","blue","green","orange")
alpha = c(.5,1,1.5,9.8)
beta = c(.5,1,1.5,1)
alpha.star = alpha + sum(horsekick)
beta.star = n + beta
plot(x.grid,x.grid,type="n",ylab="Posterior Density",
     main="Posterior Gamma Distributions", las=1, xlab = expression(lambda), ylim = c(0,1))
for(i in 1:length(alpha)){
  lines(x.grid, dgamma(x.grid, alpha.star[i], beta.star[i]), 
        type = 'l', col=colors[i],lwd=2)
}
legend("topleft", legend=c("alpha = beta = 0.5", "alpha = beta = 1", "alpha = beta = 1.5", "alpha = 9.8, beta = 1"),
       lwd=rep(2,5), col=colors, bty="n", ncol=2)
grid = seq(0.000001,0.000015,length=1000)
plot(grid, dgamma(x.grid, 9.8, 1), type = 'l')
N=10^6
alpha = 1
beta = 1


plot(grid,grid,type="n",xlim=c(0.000005,0.000015),ylim=c(0,1),ylab="Posterior Density",
     main="Posterior Distributions", las=1, xlab = expression(lambda))
alpha.star2 = sum(horsekick) + alpha
beta.star2 = 20*N + beta - sum(horsekick)

plot(grid, dbeta(grid, alpha.star2[1], beta.star2[1]), type = 'l', 
     ylab = "Posterior Density", main = "Posterior Beta Distribution", xlab = expression(theta), col = colors[1])
for(i in 1:length(alpha.star2)){
  lines(grid, dbeta(grid, alpha.star2[i], beta.star2[i]), 
        type = 'l', col=colors[i],lwd=2)
}
alpha
beta
legend("left", legend=c("alpha = beta = 0.5", "alpha = beta = 1", "alpha = beta = 20", "alpha = beta = 60"),
       lwd=rep(2,5), col=colors, bty="n", ncol=1)

beta.star2[1]

plot(grid, dbeta(grid, alpha.star2[3], beta.star2[3]), type = 'l', 
     ylab = "Posterior Density", main = "Posterior Distribution", xlab = expression(theta))
lines(grid, dbeta(grid, alpha.star2[2], beta.star2[2]), type = 'l', 
     ylab = "Posterior Density", main = "Posterior Distribution", xlab = expression(theta))
lines(grid, dbeta(grid, alpha.star2[4], beta.star2[4]), type = 'l', 
     ylab = "Posterior Density", main = "Posterior Distribution", xlab = expression(theta))


#### PART 2 ####
x = horsekick
grid = seq(0,2,.01)


alpha = c(.5,5,1,2)
beta = c(.5,1,3,2)


plot(grid,grid,type="n",xlim=c(0,1),ylim=c(0,4),xlab="",ylab="Prior Density",
     main="Prior Distributions", las=1)
for(i in 1:length(alpha)){
  prior = dbeta(grid,alpha[i],beta[i])
  lines(grid,prior,col=colors[i],lwd=2)
}

legend("topleft", legend=c("Beta(0.5,0.5)", "Beta(5,1)", "Beta(1,3)", "Beta(2,2)"),
       lwd=rep(2,5), col=colors, bty="n", ncol=3)

n = length(horsekick)
# the log-likelihood is a function of lambda and the data
poissonLogLik = function(x,theta){
  # total number of observations
  n = length(x)
  # equation
  logLik = (sum(x)*log(theta)-n*theta - sum(log(factorial(x))))
  return(logLik)
}



poi_ll = poissonLogLik(horsekick, mean(horsekick))
k = 1
BIC_M1 = -2 * poi_ll + k * log(length(horsekick))

binomLogLik = function(x, theta){
  n = length(x)
  logLik = sum(x) * log(theta) + (n*10^6 - sum(x)) * log(1 - theta) + sum(log(choose(10^6,x)))
  return(logLik)
}

bin.mle = sum(horsekick)/(length(horsekick) * 10^6)
bin_ll = binomLogLik(horsekick, bin.mle)

BIC_M2 = -2 * bin_ll + log(length(horsekick))

BIC = c(BIC_M1, BIC_M2)



set.seed(1)

S=10000
alpha = 1
beta = 1
n = length(horsekick)
post.beta.sample = rbeta(S, alpha + sum(horsekick), n*10^6+ beta - sum(horsekick))


hlp<-NULL
for(t in 1:S){
  hlp[t]<-binomLogLik(horsekick, post.beta.sample[t])
}
lph<-binomLogLik(horsekick, mean(post.beta.sample))
pdic=2*(lph-mean(hlp))
DIC_M1=-2*lph+2*pdic
alpha = 9.8
beta = 1
n = 10
post.gamma.sample = rgamma(S, alpha + sum(horsekick), n+beta)


hlp<-NULL
for(t in 1:S){
  hlp[t]<-poissonLogLik(horsekick, post.gamma.sample[t])
}
lph<-poissonLogLik(horsekick, mean(post.gamma.sample))
pdic=2*(lph-mean(hlp))
DIC_M2=-2*lph+2*pdic
DIC_M2 =121.18
DIC_M1 = 121.17
DIC = c(DIC_M1, DIC_M2)

#### Gelfand and Ghosh Binomial ####
### Calculate Gelfand and Ghosh#####################
#first create predicted values for each of our posterior samples
#n columns (as many columns as we have data points) and S rows (number of posterior samples)
b.pred_values=matrix(0,S,20)
for(i in 1:length(post.beta.sample)){
  #obtain one prediction for each hospital at each posterior sample
  b.pred_values[i,]=rbinom(rep(1,20), 10^6, rep(post.beta.sample[i], length(horsekick)))
}
#G term of gelfand and ghosh
b.g=sum((apply(b.pred_values, 2,mean)-horsekick)^2)
b.p=sum(apply(b.pred_values,2,var))
gg_criterion_bin=b.g+b.p



#### Gelfand and Ghosh Poisson ####
### Calculate Gelfand and Ghosh#####################
#first create predicted values for each of our posterior samples
post.sample.gamma = rgamma(S, 9.8 + sum(horsekick), 1 + length(horsekick))
#n columns (as many columns as we have data points) and S rows (number of posterior samples)
g.pred_values=matrix(0,S,20)
for(i in 1:length(post.sample.gamma)){
  #obtain one prediction for each hospital at each posterior sample
  g.pred_values[i,]=rpois(rep(1,20), rep(post.sample.gamma[i], length(horsekick)))
}

#G term of gelfand and ghosh
g.g=sum((apply(g.pred_values, 2,mean)-horsekick)^2)
g.p=sum(apply(g.pred_values,2,var))
gg_criterion_pois=g.g+g.p

GGH = c(gg_criterion_pois, gg_criterion_bin)

table = data.frame(BIC, DIC, GGH)
names(table) = c("BIC", "DIC", "Gelfand and Ghosh")

library(xtable)
xtable(table)
#### Bayes Factor ####

n = length(horsekick)
alpha = 9.8
beta = 1
log.m1 = alpha * log(beta) - sum(log(factorial(horsekick))) - lgamma(alpha) + 
  lgamma(sum(horsekick) + alpha) - (sum(horsekick) + alpha) * log(n + beta)

test.lm1 = lgamma(sum(horsekick) + alpha) - (sum(horsekick) + alpha) * log(n + beta)

N = 10^6
alpha = 1
beta = 1
log.m2 = sum(lchoose(N, horsekick)) + lgamma(alpha + beta) + 
  lgamma(sum(horsekick) + alpha) + lgamma(20*N + beta - sum(horsekick)) - 
  lgamma(alpha) - lgamma(beta) - lgamma(alpha + 20*N + beta)


test.lm2 = lgamma(sum(horsekick) + alpha) + lgamma(20*N + beta - sum(horsekick)) -
  lgamma(alpha) - lgamma(beta) - lgamma(alpha + 20*N + beta)

B.factor = exp(log.m1-log.m2)
exp(log.m1)/ exp(log.m2)

exp(log.m2)/exp(log.m1)


