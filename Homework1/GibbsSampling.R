### 2. Gibbs Sampling Method
dat = c(29.39, 7.94, -2.75, 6.82, -0.64, 0.63, 18.01, 12.16)
sigma_dat = c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)
pretty_table = data.frame(LETTERS[1:length(dat)], dat, sigma_dat)
names(pretty_table) = c("School" , "Est. Treatment effect (y_j)", "Est. SE (sigma_j)" )
xtable(pretty_table)

sigma_sqr_dat = sigma_dat^2
J = length(dat)
seq(A:H)
# output = NULL
# # set current values
# tau_cur = 1
# theta_cur = rep(0, length(dat))
# mu_cur = 0
# 
# 
# num_iters=8000
# # matrix/vector to store samples
# theta_samps = matrix(NA, ncol = length(dat), nrow = num_iters)
# tau_samps = rep(NA, num_iters)
# mu_samps = rep(NA, num_iters)
# # begin gibbs
# for(i in 1:num_iters){
#   tau_cur_sqr = tau_cur^2
#   theta_mean = (dat*tau_cur_sqr + mu_cur*sigma_sqr_dat)/(tau_cur_sqr + sigma_sqr_dat)
#   theta_var = tau_cur_sqr*sigma_sqr_dat/(tau_cur_sqr + sigma_sqr_dat)
#   # store next samples in chain
#   theta_samps[i,] = rnorm(J, mean = theta_mean, sd = sqrt(theta_var))
#   mu_samps[i] = rnorm(1, mean(theta_cur), sd=sqrt(tau_cur_sqr/J))
#   tau_samps[i] = sqrt(1/rgamma(1, shape = J/2-1, rate = sum((theta_cur-mu_cur)^2)/2))
# 
#   # update current values
#   theta_cur = theta_samps[i,]
#   mu_cur = mu_samps[i]
#   tau_cur = tau_samps[i]
# }
# # # discard burnin
# # mu_samps = tail(mu_samps,4000)
# # tau_samps = tail(tau_samps,4000)
# # theta_samps = tail(theta_samps,4000)
# 
# # thin
# theta_samps = theta_samps[seq(from = 1000, to = num_iters, by=2), ]
# tau_samps = tau_samps[seq(from = 1000, to = num_iters, by=2)]
# mu_samps = mu_samps[seq(from = 1000, to = num_iters, by=2)]
# 
# plot.ts(mu_samps)
# plot.ts(tau_samps)
plot.ts(theta_samps[,1])
abline(h = mean(theta_samps[,1]), col = "red")
# 


for (i in 1:num_iters) {
  # First sample theta vec
  tau_cur_sqr = tau_cur^2
  theta_mean = (dat*tau_cur_sqr + mu_cur*sigma_sqr_dat)/(tau_cur_sqr + sigma_sqr_dat)
  theta_var = tau_cur_sqr*sigma_sqr_dat/(tau_cur_sqr + sigma_sqr_dat)
  theta_cur = rnorm(J, theta_mean, sd = sqrt(theta_var))
  
  mu_cur = rnorm(1, mean(theta_cur), sd=sqrt(tau_cur_sqr/J))
  
  tau_cur = sqrt(1/rgamma(1, shape = J/2-1, rate = sum((theta_cur-mu_cur)^2)/2))
  output = rbind(output, c(theta_cur, mu_cur, tau_cur))
}
theta_cur

dim(output)
colnames(output) = c(paste("theta_", 1:J, sep=""), "mu", "tau")

output_thin = output[seq(2000, num_iters, by=2), ]
# check convergence of theta 1
plot.ts(output[,1])
abline(h = mean(output[,1]), col = "red")

summary_stats = function(x){c(mean(x), sd(x), quantile(x, c(.025,.975)))}

summ_table = matrix(unlist(apply(output_thin, 2, summary_stats)), nrow=10, byrow=T)


#### Meta Analysis
tau_opt_result = optim(1, function(x){tau_log_posterior(x, dat, sigma_sqr_dat)}, method="Brent", lower=0, upper=300, hessian=T, control=list(fnscale=-1))
tau_post_mode = tau_opt_result$value
tau_post_var = 1/(-tau_opt_result$hessian)

tau_grid = seq(1e-5, 0.5, length.out=1000)
tau_grid_log_den = sapply(tau_grid, function(x){tau_log_posterior(x, dat, sigma_sqr_dat)})
tau_post_discrete = exp(tau_grid_log_den-max(tau_grid_log_den))
sample_size = 2000
tau_sample = sample(tau_grid, sample_size, prob=tau_post_discrete, replace=T)

par(mfrow = c(1,2))
plot(tau_post_discrete~tau_grid, type="l")
hist(tau_sample)
