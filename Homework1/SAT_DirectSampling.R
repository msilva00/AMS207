# Read In Data #
library(RCurl)
# script <- getURL("https://raw.githubusercontent.com/msilva00/AMS207/master/FunctionLib/FunctionLib.R", ssl.verifypeer = FALSE)
eval(parse(text = script))


### 1. Direct Sampling
mean(dat)
dat = c(29.39, 7.94, -2.75, 6.82, -0.64, 0.63, 18.01, 12.16)
sigma_dat = c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)
sigma_sqr_dat = sigma_dat^2

tau_log_posterior = function(tau, dat, sigma_sqr_dat){
  tau_sqr = tau^2
  v_mu = 1/sum(1/(sigma_sqr_dat+tau_sqr))
  part_1 = 0.5*log(v_mu)
  part_2 = sum(-0.5*log(sigma_sqr_dat+tau_sqr))
  mu_hat = sum(dat/(sigma_sqr_dat + tau_sqr))*v_mu
  part_3 = sum(-0.5*((dat-mu_hat)^2/(sigma_sqr_dat+tau_sqr)))
  return(part_1 + part_2 + part_3)
}



tau_opt_result = optim(1, function(x){tau_log_posterior(x, dat, sigma_sqr_dat)}, method="Brent", lower=0, upper=300, hessian=T, control=list(fnscale=-1))
tau_post_mode = tau_opt_result$value
tau_post_var = 1/(-tau_opt_result$hessian)

tau_grid = seq(1e-5, 30, length.out=1000)
tau_grid_log_den = sapply(tau_grid, function(x){tau_log_posterior(x, dat, sigma_sqr_dat)})
tau_post_discrete = exp(tau_grid_log_den-max(tau_grid_log_den))
sample_size = 2000
tau_sample = sample(tau_grid, sample_size, prob=tau_post_discrete, replace=T)

par(mfrow = c(1,2))
plot(tau_post_discrete~tau_grid, type = 'l', xlab = "Grid Values", main = expression(tau ~ " Marginal Posterior Density"),
     ylab = "Density")
hist(tau_sample, main = expression(tau ~ "Marginal Posterior Histogram"), freq = F, ylab = "", xlab = "")


sample_mar_mu = function(tau, sigma_sqr_dat, dat){
  tau_sqr = tau^2
  v_mu = 1/sum(1/(sigma_sqr_dat+tau_sqr))
  mu_hat = sum(dat/(sigma_sqr_dat + tau_sqr))*v_mu
  mu = rnorm(1, mu_hat, sd = sqrt(v_mu))
  return(mu)
}

mu_sample = sapply(tau_sample, function(x){sample_mar_mu(x, sigma_sqr_dat, dat)})
hist(mu_sample, main = expression("marginal" ~ mu), freq = F, xlab = "")


sample_theta_j_post = function(mu_sample, tau_sample, dat_j, sigma_j){
  top_1 = dat_j*tau_sample^2
  top_2 = mu_sample*sigma_j
  bottom = tau_sample^2 + sigma_j
  mu = (top_1 + top_2)/bottom
  var = tau_sample^2*sigma_j/(sigma_j + tau_sample^2)
  sample_size = length(mu_sample)
  return(rnorm(sample_size, mu, sd = sqrt(var)))
}

summary_stats = function(x){c(mean(x), sd(x), quantile(x, c(.025,.975)))}

theta_sample = lapply(as.list(1:length(dat)), function(x){sample_theta_j_post(mu_sample, tau_sample, dat[x], sigma_sqr_dat[x])})

hist(theta_sample[[2]], main = "Histogram of 2000 simulated draws", freq = F, xlab = expression(theta[2]), breaks = 15)

theta_post_summ = lapply(theta_sample, summary_stats)
plot.ts(theta_sample[[1]])
abline(h = mean(theta_sample[[1]]), col = "red")
output_table = matrix(unlist(theta_post_summ), nrow=length(dat), byrow=T)

# format for latex
colnames(output_table) = c("mean", "sd", "0.025", "0.975")
row.names(output_table) =  NULL
library(xtable)
xtable(output_table)


draw_expected_theta = function(dat, sigma_sqr_dat, sample_size){
  tau_plot = seq(0.01, 30, length.out = 30)
  expected_theta_given_tau = c()
  for(tau in tau_plot){
    mu_given_this_tau = replicate(sample_size, sample_mar_mu(tau, sigma_sqr_dat, dat))
    theta_sample = lapply(as.list(1:length(dat)), function(x){sample_theta_j_post(mu_given_this_tau, tau, dat[x], sigma_sqr_dat[x])})
    expected_theta = unlist(lapply(theta_sample, mean))
    expected_theta_given_tau = rbind(expected_theta_given_tau, expected_theta)
  }
  return(list(tau = tau_plot, theta = expected_theta_given_tau))
}

result = draw_expected_theta(dat, sigma_sqr_dat, sample_size)
expected_theta = result$theta
tau_plot_grid = result$tau

par(mfrow=c(1,1))
col_palette = rainbow(length(dat))
plot(expected_theta[, 1]~tau_plot_grid, type="l", col=col_palette[1], 
     ylim=c(-5, 30), xlab="tau", main=expression("Expected " ~ theta[j] ~ "|" ~ tau), ylab = "")
for(i in 2:length(dat)){
  lines(expected_theta[, i]~tau_plot_grid, col=col_palette[i])
  
}

text(2,30, "School")
yloc = 30-2
for(i in 1:length(col_palette)){
  text(c(2,yloc), labels = LETTERS[i], col = col_palette[i])
  yloc = yloc - 2
  print(yloc)
}
points(c(2, yloc),pch = 19, col = "white", lwd = 8)
?text
