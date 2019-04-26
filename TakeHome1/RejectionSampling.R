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
sum_x = sum(x)
dat = BirdDat$RedtailedHawk
plot(BirdDat$RedtailedHawk ~ BirdDat$RouteCount)
year = BirdDat$years
ci = BirdDat$RouteCount
#### EDA ####
# Plot the raw data using the base plot function
plot(year,ci,xlab="Year", ylab = "Route Counts in California",
     main="Yearly Route Counts in California\n From 1968 to 2017",
     cex = 0.5, type = "p", col = "Blue", pch = 19)
plot(year, BirdDat$RedtailedHawk,xlab="Year", ylab = "Red Hawk Sightings",
     main="Yearly Counts of Red Hawk Sightings \n in California From 1968 to 2017",
     cex = 0.5, type = "p", col = "Blue", pch = 19)
plot(BirdDat$RedtailedHawk, ci ,xlab="Red Hawk Sightings", ylab = "Route Counts in California",
     main="Red Hawk Sightings per Route Count\n in California From 1968 to 2017",
     cex = 0.5, type = "p", col = "Blue", pch = 19)


hist(x,breaks = 10, freq = FALSE, main = "Histogram", xlab = "Red Hawk Sightings")
lines(density(x), col = "blue")

##### Latex output for Data #####
# library(xtable)
# xtable(DATATAB)
DATATAB = head(BirdDat[,1:3],10)
xtable(DATATAB)




#### log marginal posterior ####
log_marginal_posterior = function(t_v_vec, dat) {
  t = t_v_vec[1]
  exp_t = exp(t)
  v = t_v_vec[2]
  exp_v = exp(v)
  sum_dat = sum(dat)
  J = length(dat)
  part_1 = sum(sapply(dat, function(x) {
    lgamma(exp_t + x) - (lgamma(exp_t) + lgamma(x + 1))
  }))
  
  part_2 = J * exp_t * v - (J * exp_t + sum_dat) * log(exp_v + 
                                                         1) + t
  
  return(part_1 + part_2)
}



grid_len = 100
t_grid = seq(1.75, 2.6, length.out = grid_len)
v_grid = seq(-4.3, -3.3, length.out = grid_len)
joint_grid = expand.grid(t_grid, v_grid)
dim(joint_grid)


log_mat = matrix(apply(joint_grid, 1, function(x) log_marginal_posterior(x, 
                                                                         dat)), nrow = grid_len)

dim(log_mat)



library(mvtnorm)
optim_result = optim(c(2.2, -3.75), function(x) log_marginal_posterior(x, 
                                                                  dat), hessian = T, control = list(fnscale = -1))

optim_loc = optim_result$par
optim_var = solve(-optim_result$hessian)

log_density_ratio = function(t_v_vec, dat, optim_loc, optim_var, 
                             df) {
  part_1 = log_marginal_posterior(t_v_vec, dat)
  part_2 = mvtnorm::dmvt(t_v_vec, optim_loc, optim_var, df = df, 
                         log = T)
  return(part_1 - part_2)
}

sir_sample_v_t = function(t_loc, t_var, dat, df, sample_size) {
  num_iter_needed = 5 * sample_size
  proposal_sample = rmvt(num_iter_needed, delta = t_loc, sigma = t_var, 
                         df = df)
  log_den_ratio = apply(proposal_sample, 1, function(x) {
    log_density_ratio(x, dat, t_loc, t_var, df)
  })
  den_ratio = exp(log_den_ratio)
  
  resample_weights = den_ratio/sum(den_ratio)
  index_selected = sample(1:num_iter_needed, size = sample_size, 
                          prob = resample_weights, replace = T)
  output = proposal_sample[index_selected, ]
  return(output)
}

df = 4
sample_size = 2000
posterior_sample = sir_sample_v_t(optim_loc, 2 * optim_var, dat, 
                                  df, sample_size)
alpha_beta_sample = exp(posterior_sample)

library(scales)
#### Joint Sampling Importance Resampling for alpha and beta ######
contour((t_grid), (v_grid), exp(log_mat - max(log_mat)),# xlim = c(1.75, 2.6), ylim = c(-4.3, -3.3),
        xlab = expression(log(alpha)), ylab = expression(log(beta)), main = "Contour Plot")
points((posterior_sample[, 1]), (posterior_sample[, 2]), pch = 20,  cex = 0.5,
       col = alpha("red", 0.2))



# posterior for mu_i
param_posterior_param = lapply(as.list(dat), function(x) {
  cbind(alpha_beta_sample[, 1] + x, alpha_beta_sample[, 2] + 
          1)
})
theta_posterior_sample = lapply(param_posterior_param, function(x) {
  rgamma(sample_size, shape = x[, 1], rate = x[, 2])
})
sum_stats = function(x) {
  c(mean(x), sd(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
}

theta_posterior_summ = lapply(theta_posterior_sample, sum_stats)
present_table = matrix(unlist(theta_posterior_summ), nrow = length(dat), 
                       byrow = T)



colnames(present_table) = c("Mean", "Sd", "2.5%", 
                         "25%", "50%", "75%", "97.5%")
head(present_table)

alpha_summary = sum_stats(alpha_beta_sample[,1])
beta_stummary = sum_stats(alpha_beta_sample[,2])

ab_table = rbind(alpha_summary, beta_stummary, present_table[1:5,])
colnames(ab_table) = c("Mean", "Sd", "2.5%", 
                            "25%", "50%", "75%", "97.5%")
rownames(ab_table) = NULL
xtable(ab_table[,c(1,2,3,7)])

#### Traceplot and Histogram Analysis ####
# alpha
par(mfrow=c(1,2))
plot.ts(alpha_beta_sample[,1], main = expression("Traceplot for " ~ alpha), 
        xlab = "Iterations", ylab = "")
abline(h = mean(alpha_beta_sample[,1]), col = "red")
hist(alpha_beta_sample[,1], main = expression("Histogram for" ~ alpha), freq = F, breaks = 20,
     xlab = "")
lines(density(alpha_beta_sample[,1]), col = "red")

# beta
plot.ts(alpha_beta_sample[,2], main = expression("Traceplot for " ~ beta), 
        xlab = "Iterations", ylab = "")
abline(h = mean(alpha_beta_sample[,2]), col = "red")
hist(alpha_beta_sample[,2], main = expression("Histogram for" ~ beta), freq = F, breaks = 25,
     xlab = "")
lines(density(alpha_beta_sample[,2]), col = "red")

# theta_1 = lambda_1/c_1 (= mu_i in my paper)
plot.ts(theta_posterior_sample[[1]], 
        main = expression("Traceplot for " ~ mu[1]),
        xlab = "iteration", ylab = "")
abline(h = mean(theta_posterior_sample[[1]]), col = "red")
hist(theta_posterior_sample[[1]], 
        main = expression("Histogram for " ~ mu[1]),
        xlab = "", freq = F)
lines(density(theta_posterior_sample[[1]]), col = "red")

# theta_2
plot.ts(theta_posterior_sample[[2]], 
        main = expression("Traceplot for " ~ mu[2]),
        xlab = "iteration", ylab = "")
abline(h = mean(theta_posterior_sample[[2]]), col = "red")
hist(theta_posterior_sample[[2]], 
     main = expression("Histogram for " ~ mu[2]),
     xlab = "", freq = F)
lines(density(theta_posterior_sample[[2]]), col = "red")

plot.ts(theta_posterior_sample[[3]], 
        main = expression("Traceplot for " ~ mu[3]),
        xlab = "iteration", ylab = "")
abline(h = mean(theta_posterior_sample[[3]]), col = "red")
hist(theta_posterior_sample[[3]], 
     main = expression("Histogram for " ~ mu[3]),
     xlab = "", freq = F)
lines(density(theta_posterior_sample[[3]]), col = "red")


plot.ts(theta_posterior_sample[[4]], 
        main = expression("Traceplot for " ~ mu[4]),
        xlab = "iteration", ylab = "")
abline(h = mean(theta_posterior_sample[[4]]), col = "red")
hist(theta_posterior_sample[[4]], 
     main = expression("Histogram for " ~ mu[4]),
     xlab = "", freq = F)
lines(density(theta_posterior_sample[[4]]), col = "red")


# 
theta_mat = matrix(unlist(theta_posterior_sample), nrow = length(dat),
                   byrow = T)
theta_mat = t(theta_mat)
dim(theta_mat)
get_lambda = function(x){
  return(x/ci)
}

# randomly sample 200 subsamples from the of parameters 
# in order to produce posterior predictive replicates of y_i
randnums = order(sample(1:sample_size, 200, replace = F))
rep_mus = theta_mat[c(1, randnums),]
#### Posterior Predictive #####
# store replicates into a matrix
yrep = matrix(NA, nrow=201, ncol = length(x))
for(i in 1:NCOL(yrep)){
  yrep[i,] = rpois(length(x), lambda = rep_mus[i,])
}

# Estimate the probability of observing more than 450 Red Hawks 
# in a year that has a route count greater than 120
lamdas = get_lambda(rep_mus)
dim(lamdas)
yrepp = matrix(NA, nrow = dim(yrep)[1], ncol = dim(yrep)[2])
pair_counts = NA
for(i in 1:201){
  yrepp[i, ] = rpois(length(x), rep_mus[i,])
  pairs1 = cbind(na.omit(yrepp[i,]), rep_mus[i,])
  paires = data.frame(pairs1)
  pair_counts[i] = length(pairs1[paires$X1 > 450 & paires$X2>=120])
}
mean(pair_counts)/50 # probability = 0.435

# Model 2: Comparison of Actual vs Fitted Densities
par(mfrow = c(1,1))
plot(density(x), ylim = c(0,0.006),main = "Model 2: Comparison of \n Actual vs Fitted Densities",
     xlab = "Red Hawk Sightings")
for(i in 1:NCOL(yrep)){
  lines(density(yrep[i,]),col = alpha("blue", 0.1))
}
legend(x=-30, y=0.006, legend=c(expression(y[i]),  expression(paste(y[i], " replicates"))),
       col = c("black", "blue"), lty=1)

# replicates vs mean
T_means = colMeans(yrep, na.rm = TRUE)
hist(T_means,freq = F, breaks = 20, main = "", xlab = "")
abline(v = mean(x), lwd = 5, col = alpha("red", 0.5))
legend(min(T_means), 0.007,title="T = mean",fill=c("white", "red"), legend=c(expression(T(y[i])),  expression(T(y[i]^rep))))

#### GG ####
g=sum((apply(na.omit(yrep),2, mean)-x)^2)
p=sum(apply(na.omit(yrep),2, sd))
(gg_criterion=g+p) #2449.347
# col.means_mui = apply(mu.i,2, mean)
# col.means_sigi2 = apply(sig2.i, 2, mean)

par(mfrow=c(1,3))
autocorr.plot(alpha_beta_sample[,1], col=1, lwd=4, cex.axis=1.5, cex.lab=1.5, auto.layout = FALSE, 
              main=expression(alpha))
autocorr.plot(alpha_beta_sample[,2], col=1, lwd=4, cex.axis=1.5, cex.lab=1.5, auto.layout = FALSE, 
              main=expression(beta))
autocorr.plot(theta_posterior_sample[[1]], col=1, lwd=4, cex.axis=1.5, cex.lab=1.5, auto.layout = FALSE, 
              main=expression(mu[1]))






# keep = NA
# yrep2 = matrix(NA, nrow=201, ncol = length(x))
# for(i in 1:NCOL(yrep)){
#   mussss = rep_mus[i,]
#   keep = rep(0, length(mussss))
#   keep[mussss>120] = mussss[mussss > 120]
#   yrep2[i,] = rpois(length(x), lambda = keep)
# }
# yrep22 = na.omit(c(yrep2))
# num = length(yrep22[yrep22 > 450])
# den = length(yrep22)
# # proportion of sightings greater than 450
# num/den

