KMLP_PLM = function(x_rand, z, y_noise, K_max, m, B = 1000, nmc = 1000, 
                            lsig = 0.0001, usig = 10, tau = 10, kappa = 10, A = 100){
  mu_star = (2 * (1:K_max) - 1)/(2 * K_max)
  q = ncol(z)
  # Take the kernel to be the bump function
  mcmc = NULL
  mcmc$sigma = rep(NA, B + nmc)
  mcmc$K = rep(NA, B + nmc)
  mcmc$mu = matrix(NA, nrow = K_max, ncol =  B + nmc)
  mcmc$omega = array(NA, dim = c(m, K_max, B + nmc))
  mcmc$h_K = rep(NA, B + nmc)
  mcmc$beta = matrix(NA, nrow = q, ncol = B + nmc)
  mcmc$K[1] = K_max
  mcmc$beta[, 1] = rnorm(q, mean = 0, sd = 0.5)
  mcmc$sigma[1] = 1
  for (k in 1:K_max){
    mcmc$mu[k, 1] = runif(1, min = -1, max = 1)/(2 * K_max) + mu_star[k]
    for (alpha in 1:m){
      mcmc$omega[alpha, k, 1] = runif(1, min = -A, max = A)
    }
  }
  mcmc$h_K[1]= runif(1, min = lh/K_max, max = uh/K_max)
  ptm = proc.time()
  for (iter in 2:(B + nmc)){
    sigma = mcmc$sigma[iter - 1]
    K = mcmc$K[iter - 1]
    mu = mcmc$mu[1:K, iter - 1]
    omega = mcmc$omega[, 1:K, iter - 1]
    h_K = mcmc$h_K[iter - 1]
    beta = mcmc$beta[, iter - 1]
    # Update omega
    W_tmp = W(x_rand, mu, h_K)
    Design = matrix(NA, nrow = n, ncol = K * m)
    for (alpha in 1:m){
      Design[, (1 + (alpha - 1) * K) : (alpha * K)] = Design_mat_alpha(c(x_rand), mu_star, W_tmp, alpha)
    }
    V = solve(t(Design) %*% Design + diag(1/tau^2, K * m))
    eta = matrix(y_noise, nrow = n, ncol = 1) - z %*%  matrix(beta, nrow = q, ncol = 1)
    omega_mean = V %*% (t(Design) %*% eta)
    omega_new = rmvnorm(1, mean = t(omega_mean), sigma = V)
    while(max(abs(omega_new)) > A){
      omega_new = rmvnorm(1, mean = t(omega_mean), sigma = V) 
    }
    omega_new_mat = matrix(NA, nrow = m, ncol = K)
    for (alpha in 1:m){
      omega_new_mat[alpha, ] = omega_new[(1 + (alpha - 1) * K) : (alpha * K)]
    }
    # update beta
    V_beta = solve(t(z) %*% z + diag(kappa, q))
    linear_comp = matrix(y_noise, nrow = n, ncol = 1) - Design %*% matrix(omega_new, ncol = 1)
    beta_mean = V_beta %*% (t(z) %*% linear_comp)
    beta_new = rmvnorm(1, mean = t(beta_mean), sigma = V_beta)
    # Update mu
    mu_old = mu
    for (k in 1:K){
      mu_new = mu_old
      mu_new[k] = runif(1, min = -1, max = 1)/(2 * K) + mu_star[k]
      log_lik_new = log_lik_joint(eta, x_rand, mu_star, mu_new, omega_new_mat, h_K, sigma)
      log_lik_old = log_lik_joint(eta, x_rand, mu_star, mu_old, omega_new_mat, h_K, sigma)
      log_MH = log_lik_new - log_lik_old
      MH = min(1, exp(log_MH))
      if (runif(1, min = 0, max = 1) >= MH){ mu_new = mu_old }
      mu_old = mu_new
    }
    # update bandwidth h
    h_new = runif(1, min = lh/K, max = uh/K)
    log_lik_new = log_lik_joint(eta, x_rand, mu_star, mu_old, omega_new_mat, h_new, sigma)
    log_lik_old = log_lik_joint(eta, x_rand, mu_star, mu_old, omega_new_mat, h_K, sigma)
    log_MH = log_lik_new - log_lik_old
    MH = min(1, exp(log_MH))
    if (runif(1, min = 0, max = 1) < MH){ 
      mcmc$h_K[iter] = h_new
    }else{
      mcmc$h_K[iter] = h_K
    }
    # update sigma
    f = f_mean_p1(x_rand, mu_star, mu_old, omega_new_mat, h_K)
    GSSE = sum((eta - f)^2)
    sigma_inv = rgamma(1, shape = 1 + n/2, rate = 1 + GSSE/2)
    if ((1/usig > sigma_inv) || (1/lsig < sigma_inv)){
      sigma_inv = rgamma(1, shape = 1 + n/2, rate = 1 + GSSE/2)
    }
    mcmc$sigma[iter] = 1/sqrt(sigma_inv)
    K = mcmc$K[iter - 1]
    mcmc$K[iter] = mcmc$K[iter - 1]
    mcmc$mu[1:K, iter] = mu_new
    mcmc$omega[, 1:K, iter] = omega_new_mat
    mcmc$beta[, iter] = beta_new
    if(floor(iter/100) == iter/100){
      print(paste("Iteration: #", iter, sep = " "))
    }
  }
  print(proc.time() - ptm)
  return(mcmc)
}