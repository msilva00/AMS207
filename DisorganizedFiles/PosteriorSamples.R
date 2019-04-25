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



dat = BirdDat$RouteCount + BirdDat$RedtailedHawk
# plot(dat~temp$years)
grid_len = 100
t_grid = seq(-1, 4, length.out = grid_len)
v_grid = seq(-6, 3, length.out = grid_len)
joint_grid = expand.grid(t_grid, v_grid)
dim(joint_grid)

log_mat = matrix(apply(joint_grid, 1, function(x) log_marginal_posterior(x, 
                                                                         BirdDat$RouteCount)), nrow = grid_len)


dim(log_mat)


contours_norm = c(1e-08, 1e-06, 1e-04, 0.001, 0.01, seq(0.1, 
                                                        0.9, 0.2))
contour(t_grid, v_grid, exp(log_mat - max(log_mat)), levels = contours_norm, 
        ylim = c(-5, -1.8), xlim = c(1, 4))


library(mvtnorm)
optim_result = optim(c(2.6, -2.5), function(x) log_marginal_posterior(x, 
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

df = 3
sample_size = 2000
posterior_sample = sir_sample_v_t(optim_loc, optim_var, dat, 
                                  df, sample_size)

contours_norm = c(1e-08, 1e-06, 1e-04, 0.001, 0.01, seq(0.1, 
                                                        0.9, 0.2))
contour(t_grid, v_grid, exp(log_mat - max(log_mat)), 
        ylim = c(-5.5, -2), xlim = c(-1.5, 3.5))
points(posterior_sample[, 1], posterior_sample[, 2], pch = 20, 
       col = "blue")


alpha_beta_sample = exp(posterior_sample)
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
present_table

plot.ts(theta_posterior_sample[[1]])
plot(density(theta_posterior_sample[[4]]))
hist(theta_posterior_sample[[4]])
