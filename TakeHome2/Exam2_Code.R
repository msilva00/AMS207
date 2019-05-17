####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
library(RCurl)
setwd("~/private_wip/")
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
BirdDat <-read.csv(text=getURL(gitstring))
####################################################################################################################################
head(BirdDat)
x = BirdDat$RedtailedHawk
i_years = length(x)
sum_x = sum(x)
c_i = BirdDat$RouteCount


set.seed(1)

count = 0

# specify hyperparameters. In this case using the data.
(alpha_lam = 10000)
(beta_lam = alpha_lam/1000)

# alpha <- 1
# beta <- alpha/mean(y[1:40]) 
# gam <- 3
# delta <- gam/mean(y[-(1:40)])

# set up MCMC variables
N <- 20000
N.burn <- 2000
sample_save <- NULL
sample_save$N_i <- matrix(NA, nrow = N, ncol = i_years)
sample_save$theta_i <- matrix(NA, nrow = N, ncol = i_years)
sample_save$lambda_i <- matrix(NA, nrow = N, ncol = i_years)

# initialize chains
(lambda_curr <- rgamma(i_years, alpha_lam, beta_lam) )
(Ni_curr <- rpois(i_years,lambda_curr*c_i))
beta_theta = 500
Ni_curr-x
theta_curr <- rbeta(i_years, 1, beta_theta)

log_pdfN = function(lambda_i, c_i, N_i, y_i, theta_i){
  return( -lfactorial(y_i)- lfactorial(N_i - y_i) + N_i*log(lambda_i*c_i) 
          + N_i*log(1-theta_i))
}
i=1
count_accept = 0
count_reject = 0
# sampling
for(i in 1:N){
  lambda_curr <- rgamma(i_years, Ni_curr + alpha_lam, c_i + beta_lam)
  theta_curr <- rbeta(i_years, x + 1, Ni_curr - x + beta_theta)
  
  Ni_new <- rpois(i_years,lambda_curr*c_i)
  p_curr <- log_pdfN(lambda_curr, c_i, Ni_curr, x, theta_curr)
  p_new <- log_pdfN(lambda_curr, c_i, Ni_new, x, theta_curr)
  prop_dens = dpois(Ni_new,(lambda_curr*c_i))
  curr_sens = dpois(Ni_curr,(lambda_curr*c_i))
  
  # calculate acceptance probability and accept/reject acordingly
  # accpt.prob <- exp(p_new - p_curr) 
  accpt.prob = rep(NA,i_years)
  for(j in 1:length(accpt.prob)){
    accpt.prob[j] <- exp(p_new[j] - p_curr[j]) *exp(curr_sens[j] - prop_dens[j])
    q = min(1, accpt.prob[j])
    if(runif(1, 0,1) <= q){
      Ni_curr[j] = Ni_new[j]
      count_accept = count_accept + 1
    }
    
    else{
      Ni_curr[j] = Ni_curr[j]
      count_reject = count_reject + 1
    }
    
    
  }

  # save the current draws
  sample_save$theta_i[i,] <- theta_curr
  sample_save$lambda_i[i,] <- lambda_curr
  sample_save$N_i[i,] <- Ni_curr
}
count_accept/(count_accept+count_reject)
#### Ni Posterior Analysis ####

Ni_means = apply(sample_save$N_i ,2,mean)
Ni_quants = apply(sample_save$N_i ,2, quantile, probs = c(0.025, 0.095))
Ni_summary = data.frame(t(rbind(Ni_means, Ni_quants)))
head(Ni_summary)

library("plotrix")
plotCI(Ni_summary$Ni_means, F, ui=U, li=L)

random_Ni = sample_n(Ni_summary,5)
random_Ni <- random_Ni[ order(rownames(random_Ni)), ]
colnames(random_Ni) = c("Mean", "2.5", "95")
library(xtable)
xtable(random_Ni)

par(mfrow = c(1,3))
N_i1 = tail(sample_save$N_i[,1],N)
plot.ts(N_i1, ylab = "", xlab = "", main = expression("Traceplot for" ~ N[1]))
abline(h=mean(N_i1), col = "red")
hist(N_i1, ylab = "", xlab = "", main = expression("Histogram for" ~ N[1]),breaks = 20, 
     freq = F)
polygon(density(N_i1), col = alpha("pink",0.3), border = "red")
abline(v = quantile(N_i1,c(0.025,0.975)), lty=3)

acf(N_i1, main = "ACF")

par(mfrow = c(1,3))
N_i15 = tail(sample_save$N_i[,15],N-N.burn)
plot.ts(N_i15, ylab = "", xlab = "", main = expression("Traceplot for" ~ N[15]))
abline(h=mean(N_i15), col = "red")
hist(N_i15, ylab = "", xlab = "", main = expression("Histogram for" ~ N[15]),breaks = 20, 
     freq = F, ylim = c(0,max(density(N_i15)$y)))
polygon(density(N_i15), col = alpha("pink",0.3), border = "red")
abline(v = quantile(N_i15,c(0.025,0.975)), lty=3)
acf(N_i15, main = "ACF")


N_i31 = tail(sample_save$N_i[,31],N-N.burn)
plot.ts(N_i31, ylab = "", xlab = "", main = expression("Traceplot for" ~ N[31]))
abline(h=mean(N_i31), col = "red")
hist(N_i31, ylab = "", xlab = "", main = expression("Histogram for" ~ N[31]),breaks = 20, 
     freq = F, ylim = c(0,max(density(N_i31)$y)))
polygon(density(N_i31), col = alpha("pink",0.3), border = "red")
abline(v = quantile(N_i31,c(0.025,0.975)), lty=3)
acf(N_i31, main = "ACF")

par(mfrow = c(1,3))
N_i40 = tail(sample_save$N_i[,40],N-N.burn)
plot.ts(N_i40, ylab = "", xlab = "", main = expression("Traceplot for" ~ N[40]))
abline(h=mean(N_i40), col = "red")
hist(N_i40, ylab = "", xlab = "", main = expression("Histogram for" ~ N[40]),breaks = 20, 
     freq = F, ylim = c(0,max(density(N_i40)$y)))
polygon(density(N_i40), col = alpha("pink",0.3), border = "red")
abline(v = quantile(N_i40,c(0.025,0.975)), lty=3)
acf(N_i40, main = "ACF")

N_i48 = tail(sample_save$N_i[,48],N-N.burn)
plot.ts(N_i48, ylab = "", xlab = "", main = expression("Traceplot for" ~ N[48]))
abline(h=mean(N_i48), col = "red")
hist(N_i48, ylab = "", xlab = "", main = expression("Histogram for" ~ N[48]),breaks = 20, 
     freq = F, ylim = c(0,max(density(N_i48)$y)))
lines(density(N_i48), col = "red")
acf(N_i48, main = "ACF")

#### Lambdai Posterior Analysis ####

Lambdai_means = apply(sample_save$lambda_i ,2,mean)
Lambdai_quants = apply(sample_save$lambda_i ,2, quantile, probs = c(0.025, 0.095))
Lambda_summary = data.frame(t(rbind(Lambdai_means, Lambdai_quants)))
random_Lambda = sample_n(Lambda_summary,5)
random_Lambda <- random_Lambda[ order(rownames(random_Lambda)), ]
colnames(random_Lambda) = c("Mean", "2.5", "95")
xtable(random_Lambda[-3,])

par(mfrow = c(1,3))
Lambda_i17 = tail(sample_save$lambda_i[,17],N)
plot.ts(Lambda_i17, ylab = "", xlab = "", main = expression("Traceplot for" ~ lambda[17]))
abline(h=mean(Lambda_i17), col = "red")
hist(Lambda_i17, ylab = "", xlab = "", main = expression("Histogram for" ~ lambda[17]),breaks = 20, 
     freq = F)
lines(density(Lambda_i17), col = "red")
acf(Lambda_i17, main = "ACF")

par(mfrow = c(1,3))
Lambda_i26 = tail(sample_save$lambda_i[,26],N)
plot.ts(Lambda_i26, ylab = "", xlab = "", main = expression("Traceplot for" ~ lambda[26]))
abline(h=mean(Lambda_i26), col = "red")
hist(Lambda_i26, ylab = "", xlab = "", main = expression("Histogram for" ~ lambda[26]),breaks = 20, 
     freq = F, ylim = c(0,max(density(Lambda_i26)$y)))
lines(density(Lambda_i26), col = "red")
acf(Lambda_i26, main = "ACF")

par(mfrow = c(1,3))
Lambda_i45 = tail(sample_save$lambda_i[,45],N)
plot.ts(Lambda_i45, ylab = "", xlab = "", main = expression("Traceplot for" ~ theta[45]))
abline(h=mean(Lambda_i45), col = "red")
hist(Lambda_i45, ylab = "", xlab = "", main = expression("Histogram for" ~ theta[45]),breaks = 20, 
     freq = F)
lines(density(Lambda_i45), col = "red")
acf(Lambda_i45, main = "ACF")

par(mfrow = c(1,3))
Lambda_i14 = tail(sample_save$lambda_i[,14],N)
plot.ts(Lambda_i14, ylab = "", xlab = "", main = expression("Traceplot for" ~ theta[14]))
abline(h=mean(Lambda_i14), col = "red")
hist(Lambda_i14, ylab = "", xlab = "", main = expression("Histogram for" ~ theta[14]),breaks = 20, 
     freq = F)
lines(density(Lambda_i14), col = "red")
acf(Lambda_i14, main = "ACF")

#### Posterior Predictive #####
randnums = order(sample(1:N, 200, replace = F))

rep_thets = sample_save$theta_i[c(1, randnums),]
rep_lams = sample_save$lambda_i[c(1, randnums),]
rep_Ni = sample_save$N_i[c(1, randnums),]

yrep = matrix(NA, nrow=200, ncol = length(x))
for(i in 1:200){
  yrep[i,] = rep_lams[i,]*rep_thets[i,]*c_i
}

yyy = c(yrep)
length(yyy[yyy>450])
# plot(density(x), ylim = c(0,0.006),main = "Model 2: Comparison of \n Actual vs Fitted Densities",
#      xlab = "Red Hawk Sightings")

for(i in 1:NCOL(yrep)){
  lines(density(yrep[i,]),col = alpha("red", 1))
}
legend(x=-30, y=0.006, legend=c(expression(y[i]),  expression(paste(y[i], " replicates"))),
       col = c("black", "red"), lty=1)


#### GG ####
g2=sum((apply(na.omit(yrep),2, mean)-x)^2)

p2=sum(apply(na.omit(yrep),2, sd))
(gg_criterion2=g2+p2) #2449.347

Ni_means = apply(tail(sample_save$N_i,1000) ,2,mean)
Thetai_means = apply(tail(sample_save$theta_i,1000) ,2,mean)
Lambdai_means =apply(tail(sample_save$lambda_i,1000) ,2,mean)


#### DIC ####
post_sample = cbind(Ni_means, Thetai_means, Lambdai_means)

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


m2_DIC = calculateDIC(x, post_sample, log_pdfN)$DIC #didnt work



