KMLP_fast = function(x_rand, y_noise, K_max, m, B = 1000, nmc = 1000, a_sig = 2, b_sig = 2, print_iter = FALSE){
  mu_star = (2 * (1:K_max) - 1)/(2 * K_max)
  # Take the kernel to be the bump function
  mcmc = NULL
  mcmc$sigma = rep(NA, B + nmc)
  mcmc$f_pos = matrix(NA, nrow = n, ncol = nmc)
  f_pos = matrix(NA, nrow = n, ncol = B + nmc)
  mcmc$K = rep(K_max, B + nmc)
  mcmc$omega = array(NA, dim = c(m, K_max, B + nmc))
  mcmc$h_K = rep(2/K_max, B + nmc)
  # mcmc$mu[, , 1] = matrix(runif(p * mcmc$K[1], min = x_min, max = x_max), nrow = p)
  mcmc$sigma[1] = sd(y_noise)
  # mcmc$mu[, 1:mcmc$K[1], 1] = matrix(seq(x_min + 0.1, x_max - 0.1, length = mcmc$K[1]), nrow = p)
  # mcmc$theta[1:mcmc$K[1], 1] = runif(mcmc$K[1], min = -A, max = A)
  ptm = proc.time()
  K = mcmc$K[1]
  h_K = mcmc$h_K[1]
  n = length(y_noise)
  tau = n
  # Update omega
  W_tmp = W(x_rand, mu_star, h_K)
  Design = matrix(NA, nrow = n, ncol = K * m)
  for (alpha in 1:m){
    Design[, (1 + (alpha - 1) * K) : (alpha * K)] = Design_mat_alpha(c(x_rand), mu_star, W_tmp, alpha)
  }
  H_mat = Design %*% solve (t(Design) %*% Design) %*% t(Design)
  V = solve(tau^2 * Design %*% t(Design) + diag(1, n))
  Lambda = solve(t(Design) %*% Design + diag(1/tau^2, ncol(Design)))
  GSSE = matrix(y_noise, nrow = 1) %*% V %*% matrix(y_noise, ncol = 1)
  for (iter in 1:(B + nmc)){
    omega_mean = Lambda %*% t(Design) %*% matrix(y_noise, nrow = n, ncol = 1)
    omega_new_mat = omega_mean_mat = matrix(NA, nrow = m, ncol = K)
    for (alpha in 1:m){
      omega_mean_mat[alpha, ] = omega_mean[(1 + (alpha - 1) * K) : (alpha * K)]
    }
    # update sigma
    sigma_inv = rgamma(1, shape = 1 + n/2, rate = 1 + GSSE/2)
    mcmc$sigma[iter] = 1/sqrt(sigma_inv)
    # update omega
    omega_new = rmvnorm(1, mean = t(omega_mean), sigma = mcmc$sigma[iter]^2 * Lambda)
    for (alpha in 1:m){
      omega_new_mat[alpha, ] = omega_new[(1 + (alpha - 1) * K) : (alpha * K)]
    }
    mcmc$omega[, 1:K, iter] = omega_new_mat
    # mcmc$h_K[iter] = mcmc$h_K[iter - 1]
    # 
    mu = mcmc$mu[, iter]
    K = mcmc$K[iter]
    omega = matrix(mcmc$omega[, 1:K, iter], nrow = m, ncol = K)
    h_K = mcmc$h_K[iter]
    f_pos[, iter] = f_mean_p1(x, mu_star, mu_star, omega, h_K)
    if (print_iter == TRUE){
      if (floor(iter/100) == iter/100){
        print(paste("Iteration: #", iter, sep = " "))
      }
    }
  }
  mcmc$f_pos = f_pos[, (B + 1):(B + nmc)]
  # print(proc.time() - ptm)
  return(mcmc)
}