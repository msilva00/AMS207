setwd("AMS207/AMS207_Takehome1/")
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv?token=ACRTZNDHONML5AQ2X6GYYMS4YEHGI"
BirdDat <- read.csv(text=getURL(gitstring))

x = BirdDat$RedtailedHawk
n = length(x)
sum_x = sum(x)

N = 3000
sample = c()
sample$mu_i = matrix(NA,nrow = N, ncol = length(x))
sample$beta = rep(NA, length(x)) 
mu_curr = x
beta_curr = 1
alpha = 10
a = 1
b = 1


for(i in 1:N){
  mu_curr = rgamma(length(x), x + alpha, beta_curr + 1)
  beta_curr = rgamma(1, alpha*n + a, sum(mu_curr) + b)
  sample$mu_i[i,] = mu_curr
  sample$beta[i] = beta_curr
}
dim(sample$mu_i)
mus = tail(sample$mu_i, 2000)


betas = tail(sample$beta, 2000)
mu1 = mus[,1]
mu10 = mus[,10]
mu30 = mus[,30]
mu49 = mus[,49]
par(mfrow=c(1,2))
plot.ts(mu1, xlab = "Iteration", main = expression("Traceplot of " ~ mu[1]), ylab = "")
abline(h=mean(mu1), col = "red")
hist(mu1, main = expression("Histogram of " ~ mu[1]), xlab = "", freq = F, breaks = 20)
lines(density(mu1), col = "red")

plot.ts(mu10, xlab = "Iteration", main = expression("Traceplot of " ~ mu[10]), ylab = "")
abline(h=mean(mu10), col = "red")
hist(mu10, main = expression("Histogram of " ~ mu[10]), xlab = "", freq = F, breaks = 20)
lines(density(mu10), col = "red")

plot.ts(mu30, xlab = "Iteration", main = expression("Traceplot of " ~ mu[30]), ylab = "")
abline(h=mean(mu30), col = "red")
hist(mu30, main = expression("Histogram of " ~ mu[30]), xlab = "", freq = F, breaks = 20)
lines(density(mu30), col = "red")


plot.ts(mu49, xlab = "Iteration", main = expression("Traceplot of " ~ mu[49]), ylab = "")
abline(h=mean(mu49), col = "red")
hist(mu49, main = expression("Histogram of " ~ mu[49]), xlab = "", freq = F, breaks = 20)
lines(density(mu49), col = "red")

plot.ts(betas, xlab = "Iteration", main = expression("Traceplot of " ~ beta), ylab = "")
abline(h=mean(betas), col = "red")
hist(betas, main = expression("Histogram of " ~ beta), xlab = "", freq = F, breaks = 20)
lines(density(betas), col = "red")

sum_stats = function(x) {
  c(mean(x), sd(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
}

musumary = data.frame(mu1, mu10, mu30,mu49)
musumstats = lapply(musumary, sum_stats)

mu1_quants = musumstats$mu1
mu10_quants = musumstats$mu10
mu30_quants = musumstats$mu30
mu49_quants = musumstats$mu49
beta = sum_stats(betas)
quants = rbind(mu1_quants,mu10_quants,mu30_quants,mu49_quants,beta)
colnames(quants) = c("Mean", "Sd", "2.5%", 
                        "25%", "50%", "75%", "97.5%")
rownames(quants) = NULL
postquants = quants[,c(1,2,3,7)]
postquants
# xtable(postquants)

randnums = order(sample(1:N, 200, replace = F))

rep_mus = mus[c(1, randnums),]
#### Posterior Predictive #####
yrep = matrix(NA, nrow=200, ncol = length(x))
for(i in 1:NCOL(yrep)){
  yrep[i,] = rpois(length(x), lambda = rep_mus[i,])
}
plot(density(x), ylim = c(0,0.006),main = "Model 1: Comparison of \n Actual vs Fitted Densities",
     xlab = "Red Hawk Sightings")
for(i in 1:NCOL(yrep)){
  lines(density(yrep[i,]),col = alpha("red", 0.1))
}
legend(x=-30, y=0.006, legend=c(expression(y[i]),  expression(paste(y[i], " replicates"))),
       col = c("black", "red"), lty=1)


#### GG ####
g=sum((apply(na.omit(yrep),2, mean)-x)^2)
p=sum(apply(na.omit(yrep),2, sd))
(gg_criterion=g+p) #2449.347

T_means = colMeans(yrep, na.rm = TRUE)
hist(T_means,freq = F, breaks = 20, main = "", xlab = "")
abline(v = mean(x), lwd = 5, col = alpha("red", 0.5))
legend(min(T_means), 0.007,title="T = mean",fill=c("white", "red"), legend=c(expression(T(y[i])),  expression(T(y[i]^rep))))


autocorr.plot(sample$beta, col=1, lwd=4, cex.axis=1.5, cex.lab=1.5, auto.layout = FALSE, 
              main=expression(beta))
autocorr.plot(mu1, col=1, lwd=4, cex.axis=1.5, cex.lab=1.5, auto.layout = FALSE, 
              main=expression(mu[1]))

