###################################################
# This R script demonstrate how to use the provided R code to perform
# posterior inference of the kernel mixture of polynomial regression model
###################################################
# load required packages
library(mvtnorm)
library(gplots)
library(Rcpp)
###################################################
# Simulation setup for the kernel mixture of polynomials regression model
###################################################
sigma = 0.2
n = 1000
A = 50
lh = 1.2
uh = 2
x_min = 0
x_max = 1
lsig = 0.001
usig = 10
tau = 10
p = 1
# Draw sample points from the Volterra function
Rcpp::sourceCpp('kernel_weights_evaluation.cpp')
volterra_function = function(x, N = 200){
  n = length(x)
  f0 = rep(NA, n)
  i_seq = 1:N
  mu0 = i_seq^( - 3/2) * sin(i_seq)
  for (i in 1:n){
    f0[i] = sqrt(2) * sum(mu0 * cos((i_seq - 1/2) * pi * x[i]))
  }
  return(f0)
}
x_rand = runif(n, min = 0, max = 1)
y_noise = volterra_function(x_rand, N = 1000) + rnorm(n, mean = 0, sd = sigma)
x = seq(0, 1, length = n)
y_mat = volterra_function(x, N = 10000)
plot(x_rand, y_noise, xlim = c(0.1, 0.9), type = "n", col = "blue", lwd = 3, 
     cex.main = 2, cex.lab = 1.5, xlab = "x", ylab = "y")
points(x_rand, y_noise, type = "p", col = "green", lwd = 3)
points(x, y_mat, type = "l", col = "red", lwd = 3)
###################################################
# The following KMLP function included in the KMLP_function.R file is the classical 
# MCMC sampler of the kernel mixture of polynomials regression model with K selected
# by minimizing the deviance information criterion (DIC). The classical MCMC sampler
# may be time consuming. 
# A simplified version of the MCMC sampler is given at the end of the code, which is
# more efficient than the classical MCMC sampler. 
###################################################

###################################################
# Setup the MCMC for posterior inference
###################################################
B = 1000
nmc = 1000
m = 3
K_seq = 5:15 # The range of K is 5, 6, ..., 15 for model selection
source("KMLP_function.R")
mcmc_list = vector(mode = "list", length = length(K_seq))
# Model selection: Use DIC to estimate K
DIC = rep(NA, length(K_seq))
for (K_max in K_seq){
  mcmc_list[[K_max - min(K_seq) + 1]] = KMLP(x_rand, y_noise, K_max, m)
  mcmc = mcmc_list[[K_max - min(K_seq) + 1]]
  f_pos = SSE = matrix(NA, nrow = n, ncol = nmc)
  mu_star = (2 * (1:K_max) - 1)/(2 * K_max)
  D_theta = rep(NA, nmc)
  for (iter in (B + 1):(B + nmc)){
    mu = mcmc$mu[, iter]
    K = mcmc$K[iter]
    omega = matrix(mcmc$omega[, 1:K, iter], nrow = m, ncol = K)
    h_K = mcmc$h_K[iter]
    sigma_mc = mcmc$sigma[iter]
    f_pos[, iter - B] = f_mean_p1(x_rand, mu_star, mu, omega, h_K)
    # CPO[, iter - B] = dnorm(y_noise, mean = f_pos[, iter - B], sd = sigma_mc)
    D_theta[iter - B] = -2 * sum(dnorm(y_noise, mean = f_pos[, iter - B], sd = sigma_mc, log = TRUE))
  }
  DIC[K_max - min(K_seq) + 1] = var(D_theta)/2 + mean(D_theta)
}
###################################################
# choose K_max to be the K minimizing DIC
mcmc = mcmc_list[[which.min(DIC)]]
K_max = K_seq[which.min(DIC)]
###################################################
# Compute the posterior mean of f
f_pos = matrix(NA, nrow = n, ncol = nmc)
mu_star = (2 * (1:K_max) - 1)/(2 * K_max)
for (iter in (B + 1):(B + nmc)){
  mu = mcmc$mu[, iter]
  K = mcmc$K[iter]
  omega = matrix(mcmc$omega[, 1:K, iter], nrow = m, ncol = K)
  h_K = mcmc$h_K[iter]
  f_pos[, iter - B] = f_mean_p1(x, mu_star, mu, omega, h_K)
}
###################################################
# Plot the posterior estimate of the regression function
###################################################
# Compute the mean-squared error
MSE_kmlp = mean((y_mat - apply(f_pos, 1, mean))^2)
x_lim = c(min(x) + 0.05, max(x) - 0.05)
plot(x_rand, y_noise, xlim = x_lim, type = "n", col = "yellow", lwd = 3, xlab = "x", ylab = "eta", 
     main = "(a) Kernel Mixture of Polynomials")
points(x_rand, y_noise, xlim = x_lim, type = "p", col = "green", lwd = 3, pch = 1)
grid(nx = 6, ny = 6, col = "grey80", lwd = 3)
# Plot the point-wise credible intervals
polygon(c(x, rev(x)), c(apply(f_pos, 1, quantile, 0.025), rev(apply(f_pos, 1, quantile, 0.975))), 
        col = adjustcolor("grey50",alpha.f=0.5), border = F)
# Plot the point-wise posterior mean
points(x, apply(f_pos, 1, mean), type = "l", col = "blue", lwd = 3)
points(x, y_mat, type = "l", col = "red", lwd = 5, lty = "dotdash")
legend("topright", legend = c("True function", "Kernel mixture of polynomials", "Observation"), 
       lty = c("dotdash", "solid", "solid"), pch = c(NA, NA, 1),
       col = c("red", "blue", "green"), lwd = c(5, 5, 3), bty = "n")
legend("bottomleft", legend = paste("mean-squared error:", round(MSE_kmlp, 4)), bty = "n")

###################################################
# A simplified MCMC sampler for the kernel mixture of polynomials regression model
###################################################
####################################################
# In what follows we provide a simplified version of the MCMC sampler for posterior
# inference in the kernel mixture of polynomials regression model. The trick is to
# fix the kernel centers mu_k's and choose a relatively conservatively large K. This
# simplified sampler is much more efficient than the classical sampler used above,
# and the empirical performance is also satisfactory in practice. 
###################################################
source("KMLP_function_fast.R")
K_max = 10
mcmc_fast = KMLP_fast(x_rand, y_noise, K_max, m)
###################################################
# Compute the posterior mean of f
f_pos_fast = matrix(NA, nrow = n, ncol = nmc)
mu_star = (2 * (1:K_max) - 1)/(2 * K_max)
for (iter in (B + 1):(B + nmc)){
  K = K_max
  omega = matrix(mcmc_fast$omega[, 1:K, iter], nrow = m, ncol = K)
  h_K = mcmc_fast$h_K[iter]
  f_pos_fast[, iter - B] = f_mean_p1(x, mu_star, mu_star, omega, h_K)
}
###################################################
# Plot the posterior estimate of the regression function
###################################################
# Compute the mean-squared error
MSE_kmlp_fast = mean((y_mat - apply(f_pos_fast, 1, mean))^2)
x_lim = c(min(x) + 0.05, max(x) - 0.05)
plot(x_rand, y_noise, xlim = x_lim, type = "n", col = "yellow", lwd = 3, xlab = "x", ylab = "eta", 
     main = "(a) Kernel mixture of polynomials (simplified sampler)")
points(x_rand, y_noise, xlim = x_lim, type = "p", col = "green", lwd = 3, pch = 1)
grid(nx = 6, ny = 6, col = "grey80", lwd = 3)
# Plot the point-wise credible intervals
polygon(c(x, rev(x)), c(apply(f_pos_fast, 1, quantile, 0.025), rev(apply(f_pos_fast, 1, quantile, 0.975))), 
        col = adjustcolor("grey50",alpha.f=0.5), border = F)
# Plot the point-wise posterior mean
points(x, apply(f_pos_fast, 1, mean), type = "l", col = "blue", lwd = 3)
points(x, y_mat, type = "l", col = "red", lwd = 5, lty = "dotdash")
legend("topright", legend = c("True function", "Kernel mixture of polynomials", "Observation"), 
       lty = c("dotdash", "solid", "solid"), pch = c(NA, NA, 1),
       col = c("red", "blue", "green"), lwd = c(5, 5, 3), bty = "n")
legend("bottomleft", legend = paste("mean-squared error:", round(MSE_kmlp_fast, 4)), bty = "n")


