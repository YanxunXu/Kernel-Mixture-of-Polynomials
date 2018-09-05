###################################################
# This R script demonstrate how to use the provided R code to perform
# posterior inference of the kernel mixture of polynomial partial linear model
###################################################
# load required packages
library(mvtnorm)
library(gplots)
library(Rcpp)
###################################################
# Simulation setup for the kernel mixture of polynomials partial linear model
###################################################
sigma = 1
kappa = 10
n = 500
A = 100
p = 1
q = 8
beta0 =  c(1.0337655, 0.1346184, 0.2854013, 0.6675080, 
          0.6732253, 0.5293229, -0.5073311, -3.3942323)
z = matrix(runif(n * q, min = -1, max = 1), nrow = n, ncol = q)
p = 1
lh = 1.2
uh = 2
x_min = 0
x_max = 1
lsig = 0.01
usig = 10
tau = 10
# Generate simulated data from y = z * beta + eta(x) + noises
x = seq(x_min, x_max, length = n)
y_mat = 2.5 * exp(-x) * sin(10 * pi * x)
x_rand = matrix(runif(n, min = x_min, max = x_max), nrow = p, ncol = n)
eta_rand = 2.5 * exp(-x_rand) * sin(10 * pi * x_rand) + rnorm(n, mean = 0, sd = sigma)
plot(x_rand, eta_rand, type = "p", col = "green", lwd = 3, cex.main = 2, cex.lab = 1.5, xlab = "x", ylab = "y")
points(x, y_mat, type = "l", col = "red", lwd = 3)
grid(nx = 6, ny = 6, col = "grey80", lwd = 3)
y_reg = matrix(z %*% matrix(beta0, nrow = q, ncol = 1), nrow = 1, ncol = n)
y_noise = y_reg + eta_rand
###################################################
# Setup the MCMC for posterior inference
###################################################
B = 1000
nmc = 1000
m = 3
# Choose a conservative K_max = 10; 
# Otherwise one can use the model selection approach to estimate K by minimizing DIC
K_max = 10
K_seq = 5:15 # The range of K is 5, 6, ..., 15 for model selection
sourceCpp("kernel_weights_evaluation.cpp")
source("KMLP_PLM_function.R")
DIC = rep(NA, length(K_seq))
mcmc_list = vector(mode = "list", length = length(K_seq))
for (K_max in K_seq){
  mcmc_list[[K_max - min(K_seq) + 1]] = KMLP_PLM(x_rand, z, y_noise, K_max, m)
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
    beta = matrix(mcmc$beta[, iter], nrow = q, ncol = 1)
    f_pos[, iter - B] = f_mean_p1(x_rand, mu_star, mu, omega, h_K)
    D_theta[iter - B] = -2 * sum(dnorm(y_noise, mean = c(z %*% beta) + f_pos[, iter - B], sd = sigma_mc, log = TRUE))
  }
  # log_CPO[K_max - min(K_seq) + 1] = -sum(log(apply(CPO, 1, mean)))
  DIC[K_max - min(K_seq) + 1] = var(D_theta)/2 + mean(D_theta)
}
###################################################
# choose K_max to be the K minimizing DIC
mcmc = mcmc_list[[which.min(DIC)]]
K_max = K_seq[which.min(DIC) - min(K_seq) + 1]
###################################################
# Posterior mean of f
f_pos = matrix(NA, nrow = n, ncol = nmc)
for (iter in (B + 1):(B + nmc)){
  mu = mcmc$mu[, iter]
  K = mcmc$K[iter]
  omega = matrix(mcmc$omega[, 1:K, iter], nrow = m, ncol = K)
  h_K = mcmc$h_K[iter]
  f_pos[, iter - B] = f_mean_p1(x, mu_star, mu, omega, h_K)
}
###################################################
# Plot the posterior estimate of the nonlinear function eta(x)
###################################################
# Compute the mean-squared error
MSE_kmlp = mean((y_mat - apply(f_pos, 1, mean))^2)
x_lim = c(min(x) + 0.05, max(x) - 0.05)
plot(x_rand, y_noise, xlim = x_lim, ylim = c(min(y_noise) - 2, max(y_noise) + 4), 
     type = "n", col = "yellow", lwd = 3, xlab = "x", ylab = "eta", 
     main = "(a) Kernel Mixture of Polynomials")
points(x_rand, y_noise, type = "p", col = "green", lwd = 3, pch = 1)
grid(nx = 6, ny = 6, col = "grey80", lwd = 3)
# Plot the point-wise credible intervals
polygon(c(x, rev(x)), c(apply(f_pos, 1, quantile, 0.025), rev(apply(f_pos, 1, quantile, 0.975))), 
        col = adjustcolor("grey50", alpha.f = 0.5), border = F)
# Plot the point-wise posterior mean
points(x, apply(f_pos, 1, mean), type = "l", col = "blue", lwd = 3)
points(x, y_mat, type = "l", col = "red", lwd = 3, lty = "dotdash")
legend("topright", legend = c("True function", "Kernel mixture of polynomials", "Observation"), 
       lty = c("dotdash", "solid", "solid"), pch = c(NA, NA, 1),
       col = c("red", "blue", "green"), lwd = c(3, 3, 3), bty = "n")
legend("bottomleft", legend = paste("mean-squared error:", round(MSE_kmlp, 4)), bty = "n")

