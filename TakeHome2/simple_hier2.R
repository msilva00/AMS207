####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
BirdDat <- read.csv(text=getURL(gitstring))
####################################################################################################################################



yi = BirdDat$RedtailedHawk
n = length(yi)
c_i = BirdDat$RouteCount

N = 3000
sample = c()
sample$lambda_i = matrix(NA,nrow = N, ncol = length(yi))
sample$theta = rep(NA, length(yi)) 
lam_curr = yi
theta_curr = 1
alpha = 1
theta = 1

i=1
for(i in 1:N){
  lam_curr = rgamma(length(yi), yi + 1, c_i + theta_curr)
  theta_curr = rgamma(1, alpha + n, sum(lam_curr) + theta)
  sample$lambda_i[i,] = lam_curr
  sample$theta[i] = theta_curr
}
dim(sample$lam_i)
lambdas = tail(sample$lambda_i, 2000)
apply(lambdas, 2, mean)

thetas = tail(sample$theta, 2000)
lam1 = lambdas[,1]
lam10 = lambdas[,10]
lam30 = lambdas[,30]
lam49 = lambdas[,49]
par(mfrow=c(1,2))
plot.ts(lam1, xlab = "Iteration", main = expression("Traceplot of " ~ lambda[1]), ylab = "")
abline(h=mean(lam1), col = "red")
hist(lam1, main = expression("Histogram of " ~ lambda[1]), xlab = "", freq = F, breaks = 20)
lines(density(lam1), col = "red")

plot.ts(lam10, xlab = "Iteration", main = expression("Traceplot of " ~ lambda[10]), ylab = "")
abline(h=mean(lam10), col = "red")
hist(lam10, main = expression("Histogram of " ~ lambda[10]), xlab = "", freq = F, breaks = 20)
lines(density(lam10), col = "red")

plot.ts(lam30, xlab = "Iteration", main = expression("Traceplot of " ~ lambda[30]), ylab = "")
abline(h=mean(lam30), col = "red")
hist(lam30, main = expression("Histogram of " ~ lambda[30]), xlab = "", freq = F, breaks = 20)
lines(density(lam30), col = "red")


plot.ts(lam49, xlab = "Iteration", main = expression("Traceplot of " ~ lambda[49]), ylab = "")
abline(h=mean(lam49), col = "red")
hist(lam49, main = expression("Histogram of " ~ lambda[49]), xlab = "", freq = F, breaks = 20)
lines(density(lam49), col = "red")

plot.ts(thetas, xlab = "Iteration", main = expression("Traceplot of " ~ theta), ylab = "")
abline(h=mean(thetas), col = "red")
hist(thetas, main = expression("Histogram of " ~ theta), xlab = "", freq = F, breaks = 20)
lines(density(thetas), col = "red")

sum_stats = function(x) {
  c(mean(x), sd(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
}

lamsumary = data.frame(lam1, lam10, lam30,lam49)
lamsumstats = lapply(lamsumary, sum_stats)

lam1_quants = lamsumstats$lam1
lam10_quants = lamsumstats$lam10
lam30_quants = lamsumstats$lam30
lam49_quants = lamsumstats$lam49
theta = sum_stats(thetas)
quants = rbind(lam1_quants,lam10_quants,lam30_quants,lam49_quants,theta)
colnames(quants) = c("Mean", "Sd", "2.5%", 
                     "25%", "50%", "75%", "97.5%")
rownames(quants) = NULL
postquants = quants[,c(1,2,3,7)]
postquants
# xtable(postquants)


#### Posterior Predictive #####
randnums = order(sample(1:N, 200, replace = F))
rep_lams = lambdas[c(1, randnums),]

yrep1 = matrix(NA, nrow=200, ncol = length(x))
for(i in 1:200){
  yrep1[i,] = rpois(length(x), lambda = rep_lams[i,]*c_i)
}
line(density(x))
plot(density(x), ylim = c(0,0.006),main = "",
     xlab = "Red Hawk Sightings")
for(i in 1:200){
  lines(density(yrep1[i,]),col = alpha("blue", 1))
}
lines(density(x), ylim = c(0,0.006),main = "",
     xlab = "Red Hawk Sightings", lwd=2)
legend(x=-30, y=0.006, legend=c("True Obs.", "M1 Replicates", "M2 Replicates"),
       col = c("black", "blue", "red"), lty=1)

summary(yrep1)
#### GG ####
g1=sum((apply(na.omit(yrep1),2, mean)-x)^2)
p1=sum(apply(na.omit(yrep1),2, sd))
(gg_criterion1=g1+p1) #2449.347


li_means = apply(sample$lambda_i,2, mean)
thi_means = mean(sample$theta)

#### DIC ####
post_sample2 = cbind(li_means, thi_means)

logLikelihood = function(y, theta) {
  #Get the individual parameters out of theta.
  mu = theta[1]
  sigma2 = theta[2]
  
  #sum of log likelihoods = log of product of likelihoods
  sum( dnorm(y, mu, sqrt(sigma2), log=TRUE) )
}

calculateDIC = function(y, theta_post, llFun) {
  #Calculate L
  theta_hat = apply(theta_post, 2, mean)
  L = llFun(y, theta_hat)
  
  #Calculate P
  S = nrow(theta_post) #S = number of iterations
  #Add up the log likelihoods of each iteration
  llSum = 0
  for (s in 1:S) {
    theta_s = theta_post[s,]
    llSum = llSum + llFun(y, theta_s)
  }
  P = 2 * (L - (1 / S * llSum))
  
  #Calculate DIC
  DIC = -2 * (L - P)
  
  #Return the results
  list(DIC=DIC, P=P, L=L)
}


m1_DIC = calculateDIC(yi, post_sample2, logLikelihood)$DIC
