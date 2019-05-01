sum_stats = function(x){c(mean(x), sd(x), quantile(x, c(.025,.975)))}

gibbs_sampler_for_sat = function(dat, sigma_sqr_dat, num_iters){
  tau_cur = 1
  theta_cur = rep(0, length(dat))
  mu_cur = 0
  J = length(dat)
  output = c()
  
  for (i in 1:num_iters) {
    # First sample theta vec
    tau_cur_sqr = tau_cur^2
    theta_mean = (dat*tau_cur_sqr + mu_cur*sigma_sqr_dat)/(tau_cur_sqr + sigma_sqr_dat)
    theta_var = tau_cur_sqr*sigma_sqr_dat/(tau_cur_sqr + sigma_sqr_dat)
    theta_cur = rnorm(J, theta_mean, sd = sqrt(theta_var))
    
    mu_cur = rnorm(1, mean(theta_cur), sd=sqrt(tau_cur_sqr/J))
    
    tau_cur = sqrt(1/rgamma(1, shape = J/2-1, rate = sum((theta_cur-mu_cur)^2)/2))
    output = rbind(output, c(theta_cur, mu_cur, tau_cur))}
  
  colnames(output) = c(paste("theta_", 1:J, sep=""), "mu", "tau")
  return(output)
}

num_iters=8000
output= gibbs_sampler_for_sat(dat, sigma_sqr_dat, num_iters)
output_thin = output[seq(2000, num_iters, by=2), ]
summ_table = matrix(unlist(apply(output_thin, 2, sum_stats)), nrow=10, byrow=T)
# summary_thetas = summ_table[1:8,]
# names(summary_thetas) = c("mean", "sd", "0.025", "0.975")
# xtable(summary_thetas)
# summ_table
nu_theta = summ_table[9:10,]
xtable(nu_theta)
hist(output[,9], freq = F, main = expression("Histogram for Gibbs samples of " ~ mu), xlab = "")
hist(output[,10], freq = F, main = expression("Histogram for Gibbs samples of " ~ tau), xlab = "")
