####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
# library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
RTH_Dat <-read.csv(text=getURL(gitstring))[,-4]
gitstr2 = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome3/RedShouldered.csv"
RSH_dat = na.omit(read.csv(text=getURL(gitstr2)))

# ####################################################################################################################################
# head(RTH_Dat)
# RTH_Dat$years = as.factor(RTH_Dat$years)
# RTH_Dat$years = seq(1,length(RTH_Dat[,1]))
# colnames(RTH_Dat) = NULL
# # yi1 = RTH_Dat$RedtailedHawk
yi1 = RTH_Dat[,3]
yi2 = RSH_dat[,3]
# 
# head(RSH_dat)
# RSH_dat$Year = as.factor(RSH_dat$Year)
# RSH_dat$Year = seq(1,length(RSH_dat[,1]))
# colnames(RSH_dat) = NULL
# 
# RH_transp = cbind(t(RSH_dat), t(RTH_Dat))
# RH = t(RH_transp)
# library(dplyr)
# my_data = as_tibble(RH)
# my_data
# RH_defined = my_data %>% arrange(V1)
# RH = data.frame(RH_defined)
setwd("~/AMS207/")
# write.csv(RH, "RedHawkCombined")
# RH = read.csv("RedHawkCombined")[,-1]
# yi = RH$V3
# c = RH$V2
I =50
#### Gibbs ####
N = 10000
# empty vectors/matrices
alpha_post = rep(NA, N)
beta_post =  matrix(NA, I, N)
lambda_post =  matrix(NA, I, N )
gamma_post = rep(NA, N)
zi1_post = matrix(NA, I, N)
zi2_post = matrix(NA, I, N)
alpha_curr = 0.1
beta_curr = rep(1, I)
lam_curr = rep(3, I)
gamma_curr = rep(.14, I)
# zi1_curr = rep(0, I)
zi1_curr = sample(c(0,1), I, replace = T)
zi2_curr = 1-zi1_curr
# hyperpriors
a = 4
u_b = 1
v_b = 1
u_gamma = 1
v_gamma = 0.01
u_alpha = 0.1
v_alpha = 0.1

ci = RSH_dat[,2]
i=1

for (i in 1:N) {
  alpha_post[i] = rbeta(1, I*sum(zi1_curr+ zi2_curr) + u_alpha, 
                        sum(zi1_curr + zi2_curr) + v_alpha)
  beta_post[,i] = rgamma(I, rep(a + u_b,I), lam_curr + v_b)
  # lambda_post[,i] = rgamma(I, (yi1 + yi2)*(zi1_curr + zi2_curr)+a, ci*zi1_curr + lam_curr*ci*(yi1 + yi2)*zi1_curr + beta_curr )
  # gamma_post[i] = rgamma(1, sum(yi1 + yi2 +zi2_curr + zi1_curr) + u_gamma, sum(lam_curr*ci)*sum(zi2_curr+zi1_curr) + v_gamma)
  # lambda_post[,i] = rgamma(I, yi1*zi1_curr + yi2*zi2_curr+a, ci*zi1_curr + gamma_curr*ci*zi2_curr + beta_curr )
  # lambda_post[,i] = rgamma(I, a + yi1 + yi2, beta_curr + ci + gamma_curr*ci )
  # gamma_post[i] = rgamma(1, yi1*zi1_curr + yi2*zi2_curr + u_gamma, sum(lam_curr*ci)*sum(zi2_curr+zi1_curr) + v_gamma)
  lambda_post[,i] = rgamma(I, yi1*zi1_curr + yi2*zi2_curr +a,
                           beta_curr + zi1_curr + zi2_curr + gamma_curr*ci*(zi1_curr+zi2_curr))
  gamma_post[i] = rgamma(1, sum((yi1*zi1_curr) + (yi2 * zi2_curr)) + u_gamma,
                           sum(lam_curr*ci*(zi1_curr + zi2_curr)) + v_gamma)
 
  # wi1 = wi1_num/widen
  
  # wi2_num = (1-alpha_curr) * dpois(yi2, lam_curr*ci*gamma_curr)
  # widen = alpha_curr * dpois(yi1, lam_curr*ci) + (1-alpha_curr)*dpois(yi2, lam_curr*gamma_curr*ci)
  # wi2 = round(wi2_num/widen,1)
  
  
  zi2_post[,i] = 0
  zi1_post[,i] = 1
  
  alpha_curr = alpha_post[i]
  beta_curr = beta_post[,i]
  lam_curr = lambda_post[,i]
  gamma_curr = gamma_post[i]
  zi1_curr = zi1_post[,i]
  # zi2_curr = 1 - zi1_post[,i]
  # zi1_curr = rep(1,I)
  zi2_curr = 1 - zi1_curr
}

plot.ts(tail(alpha_post, 5000))
plot.ts(tail(beta_post[1,],1000))
plot.ts(tail(beta_post[2,],1000))
plot.ts(tail(beta_post[20,],1000))
plot.ts(tail(gamma_post,5000))
plot.ts(tail(lambda_post[1,], 5000))
plot.ts(tail(lambda_post[10,], 5000))
plot.ts(tail(lambda_post[50,], 5000))
acf(tail(lambda_post[2,], 1000))
plot.ts(zi1_post[3,])
median(tail(lambda_post[4,], 5000))

alpha = tail(alpha_post, 5000)
beta1_mean = mean(tail(beta_post[1,],1000))
mean(alpha)



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
