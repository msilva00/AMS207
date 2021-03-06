#################################################### Preliminary Data Gathering ####################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
BirdDat <- read.csv(text=getURL(gitstring))
####################################################################################################################################

x = BirdDat$RedtailedHawk
ndat = length(x)
# alpha1 = mean(x)
# beta1 = 1
mean(x)/var(x)
mean(x)^2/var(x)
(alpha1 = (mean(x))^2/var(x))

beta1 = 0.01
alpha1 = 0.01

x.grid = seq(365,385, length = 1000)
plot(x.grid, dgamma(x.grid, sum(x) + alpha1, ndat + beta1), type = "l", lty = 4, lwd=2, col = "black")
alpha = c(10,0.001,alpha1, mean(x))
beta = c(1,0.001,beta1, 1)
alpha.star = alpha + sum(x)
beta.star = length(x) + beta
xgrid1 = seq(340,385, length = 1000)
colors = c("blue","red","black", "yellow")
linetype = c(6,2,5,3)
plot(xgrid1,xgrid1,type="n",ylab="Posterior Density",
     main="Posterior Gamma Distributions", las=1, xlab = expression(lambda), ylim = c(0,0.20))
for(i in 1:length(alpha)){
  lines(xgrid1, dgamma(xgrid1, alpha.star[i], beta.star[i]), 
        type = 'l', col=colors[i],lwd=2, lty = linetype[i])
}
legend("topleft", legend=c(paste("alpha =", alpha[1], "; beta =", beta[1]), 
                           paste("alpha = beta =", alpha[2]), 
                           paste("alpha = ", round(alpha1,2), "\n beta =", round(beta1,2), "\n"),
                           paste("alpha = ", round(mean(x),1), "\n beta = 1")),
       lwd=rep(2,5), col=colors, lty = linetype, bty="n", ncol=1)
text(376,0.16, "Values almost identical")
####################################################################################################################################
#### MODEL 2 - Beta-Binomial ####

sumx = sum(x)
N = 10^6
alpha0 = 0.01
beta0 = 0.01
xgrid = seq(0.000345,0.0004,length = 1000)
plot(xgrid, dbeta(xgrid, sumx + alpha0, ndat*N + beta0 - sumx  ),
     type = "l", lty = 4, lwd=2, col = "black")


alpha = c(60, mean(x)/N,1.019, 0.01)
beta = c(60,1,1, 0.01)

alpha.star =  sumx + alpha
beta.star = ndat*N + beta - sumx 

library(scales)
colors = c("red","blue","yellow", "black")
linetype = c(5,1,4,3)
lwidth = c(2,4,3,2)
# opacity = c(1,)
# par(mar = c(2, 8, 2, 2))
plot(xgrid,xgrid,type="n",ylab="Posterior Density",
     main="Posterior Beta Distributions", las=1, xlab = expression(theta), ylim = c(0,15e04),
     xlim = c(0.000345,0.0003875), yaxt = "none")
axis(2, xaxp = c(0,signif(15e04,0),4))
for(i in 1:length(alpha)){
  lines(xgrid, dbeta(xgrid, alpha.star[i], beta.star[i]), 
        type = 'l', col=alpha(colors[i], 1),lwd=lwidth[i], lty = linetype[i])
}
# legend("topleft", legend=c("alpha = beta = 60", "alpha = 376.1, beta = 1", "alpha = beta = 1", "alpha = 0.01, beta = 0.01"),
#        lwd=rep(2,5), col=colors, lty = linetype, bty="n", ncol=1)
legend("topleft", legend=c(paste("alpha = beta = ", alpha[1]), 
                           paste("alpha =", "3.761e-04", ", beta =", beta[2]), 
                           "alpha = beta = 1", 
                           paste("alpha =", alpha[4], ", beta =", beta[4])),
       lwd=lwidth, col=colors, lty = linetype, bty="n", ncol=1)

for(i in 1:4){
  zzz = rbeta(10000,alpha[i],beta[i])
  print(paste(round(mean(zzz), 5), round(var(zzz), 5)))
}

var(x)/N
####### BIC ######
poissonLogLik = function(x,theta){
  # total number of observations
  n = length(x)
  # equation
  logLik = (sum(x)*log(theta)-n*theta - sum(lfactorial(x) ))
  return(logLik)
}
N = 10^6
binomLogLik = function(x, theta, N){
  n = length(x)
  logLik = sum(x) * log(theta) + (n*N - sum(x)) * log(1 - theta) + sum(lchoose(N, x))
  return(logLik)
}

pois_ll = poissonLogLik(x, mean(x))

BIC_M1 = -2 * pois_ll + log(ndat)

bin.mle = sum(x)/(length(x) * N)
bin_ll = binomLogLik(x, bin.mle, N)

BIC_M2 = -2 * bin_ll + log(length(x))
print(c(BIC_M1, BIC_M2))
(BIC = c(BIC_M1, BIC_M2, abs(diff(c(BIC_M1, BIC_M2)))))

#### Bayes Factor #####
n = length(x)
alpha = 0.01
beta = 0.01
log.m1 = alpha* log(beta) + lgamma(sum(x) + alpha) - lgamma(alpha) -
  (sum(x) + alpha)*log(beta + n) - sum(lfactorial(x))


N = 10^6
alpha = 3.761e-4
beta = 1
log.m2 = sum(lchoose(N,x)) + lbeta(sum(x) + alpha, n*N - sum(x) + beta) - lbeta(alpha, beta)

(B.factor = exp(log.m1-log.m2)) # close to 1 means no evidence to favor either model over the other

#### DIC Pois-Gamma ####
S=10000
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
alpha = 3.761e-04
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

(DIC = c(DIC_M1, DIC_M2, abs(diff(c(DIC_M1, DIC_M2)))))

### Using Replicates:
# #### Gelfand and Ghosh Poisson ####
# ### Calculate Gelfand and Ghosh#####################
# #first create predicted values for each of our posterior samples
# #n columns (as many columns as we have data points) and S rows (number of posterior samples)
# g.pred_values=matrix(0,S,n)
# for(i in 1:length(post.gamma.sample)){
#   #obtain one prediction for each hospital at each posterior sample
#   g.pred_values[i,]=rpois(rep(1,n), rep(post.gamma.sample[i], length(x)))
# }
# 
# #G term of gelfand and ghosh
# g.g=sum((apply(g.pred_values, 2,mean)-x)^2)
# g.p=sum(apply(g.pred_values,2,var))
# gg_criterion_pois=g.g+g.p
# 
# 
# #### Gelfand and Ghosh Binomial ####
# ### Calculate Gelfand and Ghosh#####################
# #first create predicted values for each of our posterior samples
# #n columns (as many columns as we have data points) and S rows (number of posterior samples)
# b.pred_values=matrix(0,S,n)
# for(i in 1:length(post.beta.sample)){
#   #obtain one prediction for each hospital at each posterior sample
#   b.pred_values[i,]=rbinom(rep(1,n), N, rep(post.beta.sample[i], length(x)))
# }
# 
# #G term of gelfand and ghosh
# b.g=sum((apply(b.pred_values, 2,mean)-x)^2)
# b.p=sum(apply(b.pred_values,2,var))
# gg_criterion_bin=b.g+b.p
# 
# (GGH = c(gg_criterion_pois, gg_criterion_bin, abs(diff(c(gg_criterion_pois, gg_criterion_bin)))))
# abs(diff(GGH))


#### Gelfand and Ghosh (closed form) ####
alpha = 0.01
beta = 0.01
Epois = (sum(x) + alpha)/(n+beta)
VarPois = (sum(x) + alpha)*(n + beta + 1)/(n + beta)^2
G = sum((Epois - x)^2)
P = VarPois
(GG11 = G + P)

alpha = 1.01
beta = 1
alpha_star = sum(x) + alpha
beta_star = n*N + beta - sum(x)
Ebetabin = N*alpha_star/(alpha_star + beta_star)
varnum = N*alpha_star*beta_star*(alpha_star + beta_star + n)
varden = (alpha_star + beta_star)^2*(alpha_star + beta_star + 1)
VarBetabin = varnum/varden
G2 = sum((Ebetabin - x)^2)
P2 = VarBetabin
GG2 = G2 + P2

(GGH = c(GG11, GG2,abs(diff(c(GG11,GG2))) ))








BF = c(round(B.factor,2), "        ", "        ")
(table = data.frame(BIC, DIC, GGH, BF))
names(table) = c("BIC", "DIC", "Gelfand and Ghosh", "Bayes Factor")
rownames(table) = c("M1", "M2", "Difference")
table
library(xtable)
xtable(table)


