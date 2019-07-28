#### THE REAL DEAL FINAL DRAFT SERIOUSLY DONT DELETE ME ####
####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
RTH_Dat <-read.csv(text=getURL(gitstring))[,-4]
gitstr2 = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome3/RedShouldered.csv"
RSH_dat = na.omit(read.csv(text=getURL(gitstr2)))
# 
# ####################################################################################################################################
head(RTH_Dat)
RTH_Dat$years = as.factor(RTH_Dat$years)
RTH_Dat$years = seq(1,length(RTH_Dat[,1]))
colnames(RTH_Dat) = NULL
# yi1 = RTH_Dat$RedtailedHawk
yi1 = RTH_Dat[,3]
yi2 = RSH_dat[,3]


par(mfrow = c(1,2))
hist(yi1,breaks = 10, freq = FALSE, main = "", xlab = "Red-Tailed Hawk")
lines(density(yi1), col = "blue")
hist(yi2,breaks = 10, freq = FALSE, main = "", xlab = "Red-Shoudered Hawk")
lines(density(yi2), col = "red")
title("Histograms of Red Hawk Sighting \n by Species", outer = T, line = -3)

y = c(yi1, yi2)
hist(y, breaks = 20, freq = F, main = "Histogram of Both Species of Red Hawks", ylab = "Density", xlab = "")
lines(density(y), col = "red")
plot(density(y), col = "red")
head(RSH_dat)
RSH_dat$Year = as.factor(RSH_dat$Year)
RSH_dat$Year = seq(1,length(RSH_dat[,1]))
colnames(RSH_dat) = NULL

# RH_transp = cbind(t(RSH_dat), t(RTH_Dat))
# RH = t(RH_transp)
# library(dplyr)
# my_data = as_tibble(RH)
# my_data
# RH_defined = my_data %>% arrange(V1)
# RH = data.frame(RH_defined)
# setwd("~/AMS207")
# write.csv(RH, "RedHawkCombined")
# RH = read.csv("RedHawkCombined")[,-1]
I =50
#### Gibbs ####
N = 15000
alpha_post = rep(NA, N)
beta_post =  matrix(NA, I, N)
lambda_post =  matrix(NA, I, N )
gamma_post = rep(NA, N)
zi1_post = matrix(NA, I, N)
zi2_post = matrix(NA, I, N)
alpha_curr = 0.1
beta_curr = rep(0.01, I)
lam_curr = rep(3, I)
gamma_curr = rep(0.1, I)
# zi1_curr = rep(0, I)
zi1_curr = sample(c(0,1), I, replace = T)
zi2_curr = 1-zi1_curr
a = 1
u_alpha = 100
v_alpha = 0.01
u_b = 5
v_b = 15
u_gamma = 1
v_gamma = 0.9
ci = RSH_dat[,2]
prob_ratio = rep(NA,N)
i=1
for (i in 1:N) {
  alpha_post[i] = rbeta(1, I*sum(zi1_curr+ zi2_curr) + u_alpha, sum(zi1_curr + zi2_curr) + v_alpha)
  beta_post[,i] = rgamma(I, rep(a + u_b,50), lam_curr + v_b)
  lambda_post[,i] = rgamma(I, yi1*zi1_curr + yi2*zi2_curr +a,
                           beta_curr + zi1_curr + zi2_curr + gamma_curr*ci*(zi1_curr+zi2_curr))
  gamma_post[i] = rgamma(1, sum((yi1*zi1_curr) + (yi2 * zi2_curr)) + u_gamma,
                         sum(lam_curr*ci*(zi1_curr + zi2_curr)) + v_gamma)*0.1346154
  wi1_num = alpha_curr * dpois(yi1, lam_curr*ci)
  wi1_den = wi1_num + (1-alpha_curr)*dpois(yi2, lam_curr*gamma_curr*ci)
  wi1 = wi1_num/wi1_den
  zi1_post[,i] = rbinom(I, 1, p = wi1)
  zi2_post[,i] = 1 - zi1_post[,i]
  
  alpha_curr = alpha_post[i]
  beta_curr = beta_post[,i]
  lambda_curr = lambda_post[,i]
  gamma_curr = gamma_post[i]/0.1346154
  zi1_curr = zi1_post[,i]
  zi2_curr = 1 - zi1_post[,i]
}
length(bth)
beta.interval = apply(beta_post[,1000:1049], 2, quantile, c(0.025, 0.975))
beta.mean = apply(beta_post[,1000:1049], 2,mean)
ind = 1:50
dim(beta.interval)
dim(rbind(ind,ind))
matplot(rbind(ind,ind), beta.interval, type = "l", lty = 1, col = 1, xlab = "Observed", ylab = "Temperature")
points(ind, beta.mean, pch = 20, cex = 0.8)

n.burn <- 5000
n.thin <- 2
bth <- seq(n.burn+1, N, by=n.thin)
library(scales)
plot.ts(tail(alpha_post, 5000))
plot.ts(tail(beta_post[1,],1000))
plot.ts(tail(beta_post[2,],1000))
plot.ts(tail(beta_post[20,],1000))
plot.ts(tail(gamma_post,5000))
plot.ts(tail(lambda_post[1,], 5000))
plot.ts(tail(lambda_post[10,], 5000))
plot.ts(tail(lambda_post[50,], 5000))

##### hypers spec #####
# alpha hist
plot.ts(tail(alpha_post, 5000), xlab = "iteration", ylab = expression(alpha))
plot.ts(gamma_post[bth], xlab = "iteration", ylab = expression(gamma))
par(mfrow = c(1,2))
hist(tail(alpha_post, 5000), main = expression("Histogram of "~ alpha), 
     ylab = "Density", freq = F, xlab = "")
lines(density(tail(alpha_post, 5000)), col = "red")

polygon(density(tail(alpha_post, 5000)), col = alpha("pink",0.3), border = "red")
abline(v = quantile(tail(alpha_post, 5000),c(0.025,0.975)), lty=3)

apply(zi1_post,1, mean)

# gamma hist
hist(gamma_post[bth], main = expression("Histogram of "~ gamma), 
     ylab = "Density", freq = F, xlab = "")
lines(density(gamma_post[bth]), col = "red")

polygon(density(gamma_post[bth]), col = alpha("pink",0.3), border = "red")
abline(v = quantile(gamma_post[bth],c(0.025,0.975)), lty=3)

par(mfrow = c(1,2))
acf(tail(alpha_post, 5000), main = expression("ACF for "~alpha ))
acf(gamma_post[bth],main = expression("ACF for "~gamma ))

# beta trace
par(mfrow = c(3,1)) 
par(mar = c(3,5,1,3))
plot.ts(beta_post[1,bth], ylab = expression(b[1]))
plot.ts(beta_post[30,bth],ylab = expression(b[30]))
plot.ts(beta_post[50,bth], ylab = expression(b[50]))

# beta hist
par(mfrow = c(1,2)) 
hist(beta_post[1,bth], main = expression("Histogram of "~ b[1]), 
     ylab = "Density", freq = F, xlab = "", ylim = c(0,3.5))
lines(density(beta_post[1,bth]), col = "red")

polygon(density(beta_post[1,bth]), col = alpha("pink",0.3), border = "red")
abline(v = quantile(beta_post[1,bth],c(0.025,0.975)), lty=3)

acf(beta_post[1,bth],main = expression("ACF for "~b[1] ))

plot.ts(gamma_post[bth])

# lambda trace
par(mfrow = c(3,1)) 
par(mar = c(3,5,1,3))
plot.ts(lambda_post[1,bth], ylab = expression(lambda[1]))
plot.ts(lambda_post[30,bth],ylab = expression(lambda[30]))
plot.ts(lambda_post[50,bth], ylab = expression(lambda[50]))

# lambda hist
par(mfrow = c(1,2)) 
hist(lambda_post[1,bth], main = expression("Histogram of "~ lambda[1]), 
     ylab = "Density", freq = F, xlab = "", ylim = c(0,1.3))
lines(density(lambda_post[1,bth]), col = "red")

polygon(density(lambda_post[1,bth]), col = alpha("pink",0.3), border = "red")
abline(v = quantile(lambda_post[1,bth],c(0.025,0.975)), lty=3)

acf(lambda_post[1,bth],main = expression("ACF for "~lambdab[1] ))

plot.ts(gamma_post[bth])

plot.ts(lambda_post[1,bth])
plot.ts(lambda_post[30,bth])
plot.ts(lambda_post[50,bth])
plot.ts(zi1_post[3,])
acf(tail(gamma_post,5000))
alpha = tail(alpha_post, 5000)
beta1_mean = mean(tail(beta_post[1,],1000))
mean(alpha)

library(mcmcplots)
params = rep(NA, I)
for(i in 1:I){
  params[i] = paste("bi", i, sep = "")
}



# 
# i1_years = length(yi1)
# c_i2 = RTH_Dat$RouteCount
# colnames(RTH_Dat) = NULL
# head(RSH_dat)
# yi2 = RSH_dat$RedShouldered
# i2_years = RSH_dat$Year
# c_i2 = RSH_dat$RouteCount
# colnames(RSH_dat) = NULL
# head(RSH_dat)
# 
# RH_dat = rbind(RTH_Dat[,2:3], RSH_dat[,])

#### Posterior Predictive #####
psamps = 50
randnums = sample(x=seq(5000,5000+psamps), psamps, replace = F)
rep_lams = lambda_post[,c(randnums)]
rep_bis = beta_post[,c(randnums)]
rep_alphas = alpha_post[c(randnums)]
rep_gams = gamma_post[c(randnums)]
i=1
yi1_rep = matrix(NA, nrow = psamps, ncol = 50)
yi2_rep = matrix(NA, nrow = psamps, ncol = 50)

c=0
for(i in 1:psamps){
  # for (j in 1:2) {
  #   yrep[i,j] = rep_alphas[i]*rpois(1, lambda = rep_lams[,i]*ci)  +
  #     (1-rep_alphas[i])*rpois(1, lambda = rep_gams[i]*rep_lams[,i]*ci)  
  # }
  
  yi1_rep[i,] = rep_alphas[i]*rpois(length(yi1), lambda = rep_lams[,i]*ci)+(1-rep_alphas[i])*rpois(50, lambda = rep_lams[,i]*rep_gams[i]*ci)
  yi2_rep[i,] = yi2 * (1-rep_alphas[i])*rpois(50, lambda = rep_lams[,i]*rep_gams[i]*ci)
}

yi2
plot(density(yi1))
for(i in 1:psamps){
  lines(density(yi1_rep[i,]), col = alpha("red", 0.1))
}

plot(density(yi1))
lines(density(yi1_rep[1,]), col = alpha("red", 1))
plot(density(yi2))
lines(density(yi2_rep[1,]), col = alpha("red", 1))


yrep = matrix(NA, nrow = psamps, ncol = 100)
for (i in 1:psamps) {
  yrep[i,] = c(yi1_rep[i,],yi2_rep[i,])
}

# y.interval = apply(yrep, 2, quantile, c(0.025, 0.975))
# beta.mean = apply(beta_post[,1000:1049], 2,mean)
# ind = 1:50
# dim(y.interval)
# dim(rbind(ind,ind))
# matplot(rbind(ind,ind), beta.interval, type = "l", lty = 1, col = 1, xlab = "Observed", ylab = "Temperature")
# points(ind, beta.mean, pch = 20, cex = 0.8)


plot(density(y),main = "New Observation", xlab="")
for (i in 1:psamps) {
  lines(density(yrep[i,]), col = alpha("red", .09))
}
abline(v=mean(yi2))
abline(v=mean(yi1))
density(y)
which(samples$tau == max(samples$tau))
par(mfrow=c(1,2))
plot(yi1,pch = 19, cex = 0.5, col = "blue", ylab = "Red-Tailed Counts", xlab = "year (1968-2017)")
points(apply(yi1_rep,2, mean),pch = 19, cex = 0.5, col = "green")
legend("bottomright", legend=c(expression(y[i1]),  expression(bar(y[it]))),
       col = c("blue", "green"), lty=1)

mean_yi2 = apply(yi2_rep,2, mean)
mean_yi1 = apply(yi1_rep,2, mean)
pRT = yi1/(yi1+yi2)
pRS = yi2/(yi1+yi2)
pRTc = 1-pRT
pRSc = 1-pRS

yi2_rep[50,]
plot(yi1_rep[50,]/sum(yi2_rep[50,] + yi1_rep[50,]),pch = 19, cex = 0.5, col = "darkred", 
     ylab = "", xlab = "", main = "Probability of observing RTH")
points(yi1/sum(yi1 + yi2),pch = 19, cex = 0.5, col = "yellow")
legend("bottomleft", legend=c("Predicted",  "Observed"),
       col = c("darkred", "pink"), lty=1)
