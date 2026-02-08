# ------------------------------------------------------------------------------
# Function: Model Data Generation
# ------------------------------------------------------------------------------

generate_data <- function(n, dB = 2, dM = 3, dS = 1,
                          beta_tB, beta_MB, beta_Mt,
                          beta_rB, beta_rM,
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun, pi2_fun) {
  
  X_B <- eps_B_fun(n)
  eps_tr <- eps_tr_fun(n)
  X_t <- X_B %*% beta_tB + eps_tr[, 1]
  eps_M <- eps_M_fun(n)
  X_M <-  X_B %*% t(beta_MB) + X_t %*% beta_Mt + eps_M
  X_r <- X_B %*% beta_rB + X_M %*% beta_rM + eps_tr[, 2]
  X_S <- matrix(rnorm(n * dS), n, dS)
  
  data <- data.frame(C = rep(NA, n), X_r = X_r)
  data$X_t <- X_t
  data$X_B <- I(X_B)
  data$X_M <- I(X_M)
  data$X_S <- I(X_S)
  
  # coarsen data
  pi1 <- pi1_fun(data[, c("X_B", "X_t")])
  pi2 <- pi2_fun(data[, c("X_B", "X_t", "X_M")])
  
  C <- rep(NA, n)
  for (i in 1:n) {
    p1 <- 1 - pi1[i]
    p2 <- pi1[i] * (1 - pi2[i])
    p3 <- pi1[i] * pi2[i]
    C[i] <- sample(c(1, 2, Inf), 1, prob = c(p1, p2, p3))
  }
  data$X_M[C == 1, ] <- NA
  data$X_r[C %in% c(1,2)] <- NA
  data$X_S[C %in% c(1,2)] <- NA
  data$C <- C
  
  return(data)
}


# ------------------------------------------------------------------------------
# Function: Parameter Estimation Functions
# ------------------------------------------------------------------------------

estimate_tB <- function(data) {
  solve(var(data$X_B) %*% (t(data$X_B) %*% data$X_B) , 
        var(data$X_B) %*%t(data$X_B) %*% data$X_t)
}

estimate_M <- function(data, beta_tB, pi1_fun) { 
  eps_t <- data$X_t - as.matrix(data$X_B) %*% beta_tB
  W <- ifelse(data$C >= 2, 1 / pi1_fun(data[, c("X_B", "X_t")]), 0) 
  
  fit <- lm(as.matrix(data$X_M) ~ data$X_t + as.matrix(data$X_B),
            weights = W)
  if (!is.matrix(coef(fit))) {  return(list(beta_Mt = as.matrix(coef(fit)[2]), beta_MB = t(coef(fit)[-c(1,2)]))) }
  list(beta_Mt = coef(fit)[2, ], beta_MB = t(coef(fit)[-c(1,2), ]))
}

estimate_r <- function(data, beta_tB, beta_Mt, beta_MB, pi1_fun, pi2_fun) {
  eps_t <- data$X_t - as.matrix(data$X_B) %*% beta_tB
  eps_M <- (as.matrix(data$X_M) - data$X_t %*% beta_Mt - 
              as.matrix(data$X_B) %*% t(beta_MB))
  W <- ifelse(is.infinite(data$C),
              1 / (pi1_fun(data[, c("X_B", "X_t")]) *
                     pi2_fun(data[, c("X_B", "X_t", "X_M")])),
              0)
  
  Z2 <- cbind(as.matrix(data$X_B), eps_t, data$X_M)
  #X_B_adj <- t(solve(var(data$X_B), t(as.matrix(data$X_B)))) - 
  #  t(solve(cov.wt(eps_M[data$C == Inf, ], wt = W[data$C == Inf])$cov, t(eps_M))) %*% 
  #  (beta_Mt %*% t(beta_tB) + beta_MB)
  #Z1 <- cbind(X_B_adj, eps_t, eps_M)
  Z1 <- cbind(as.matrix(data$X_B), eps_t, eps_M)
  Xr <- data$X_r
  
  # weighted cross-products
  WZ2   <- sqrt(W) * Z2
  WZ1   <- sqrt(W) * Z1
  WXr  <- sqrt(W) * Xr
  
  # solve normal equations: (Σ W_i Zl_i Z_i^T)^(-1) (Σ W_i Zl_i Xr_i)
  crossprod_mat <- crossprod(WZ1[is.infinite(data$C), ], WZ2[is.infinite(data$C), ])
  rhs_vec       <- crossprod(WZ1[is.infinite(data$C), ], WXr[is.infinite(data$C)])

  theta_hat <- solve(crossprod_mat, rhs_vec)
  return(theta_hat)
}


estimate_parameters <- function(data, pi1_fun, pi2_fun) {
  beta_tB <- estimate_tB(data)
  M_est <- estimate_M(data, beta_tB, pi1_fun)
  r_est <- estimate_r(data, beta_tB, M_est$beta_Mt, M_est$beta_MB, 
                      pi1_fun, pi2_fun)
  
  list(beta_tB = beta_tB,
       beta_Mt = M_est$beta_Mt,
       beta_MB = M_est$beta_MB,
       beta_rB = as.matrix(r_est[1:length(beta_tB)]),
       gamma   = r_est[length(beta_tB) + 1],
       beta_rM = as.matrix(r_est[-(1:(length(beta_tB) + 1))]))
}

# ------------------------------------------------------------------------------
# Checks:
#   Non-Gaussian Error  : Positive
#   Multivariate B and M: Positive
#   Correct parameters  : Positive
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: compute residuals and moments needed for EIF variance and optimal pi
# ------------------------------------------------------------------------------

compute_moments <- function(data, estimates, pi1_0_fun, pi2_0_fun) {
  X_B <- data$X_B
  X_t <- data$X_t
  X_M <- data$X_M
  X_r <- data$X_r
  C <- data$C
  beta_tB <- estimates$beta_tB
  beta_MB <- estimates$beta_MB
  beta_Mt <- estimates$beta_Mt
  beta_rB <- estimates$beta_rB
  beta_rM <- estimates$beta_rM
  gamma <- estimates$gamma
  eps_t <- X_t - X_B %*% beta_tB
  eps_M <- X_M - X_t %*% beta_Mt - X_B %*% t(beta_MB)
  eps_r_cond_t <- X_r - X_B %*% beta_rB - X_M %*% beta_rM - eps_t %*% gamma
  
  # Moments
  var_eps_t <- var(eps_t)
  Weight_2inf <- ifelse(data$C >= 2, 1 / pi1_0_fun(cbind(X_B, X_t)), 0)
  var_eps_M <- cov.wt(eps_M[data$C >= 2, , drop = FALSE], wt = Weight_2inf[data$C >= 2])$cov
  Weight_inf <- ifelse(is.infinite(C),
                       1 / (pi1_0_fun(cbind(X_B, X_t)) *
                              pi2_0_fun(cbind(X_B, X_t, X_M))), 0)
  var_eps_r_cond_t <- cov.wt(as.matrix(eps_r_cond_t[is.infinite(C)]), 
                             wt = Weight_inf[is.infinite(C)])$cov
  
  avar_const_1 <- as.numeric(t(beta_rM) %*% var_eps_M %*% beta_rM / var_eps_t^2)
  avar_const_2 <- solve(var_eps_M, beta_Mt) * c(sqrt(var_eps_r_cond_t))
  
  list( 
    eps_t = eps_t,
    eps_M = eps_M,
    Weight_2inf = Weight_2inf,
    avar_const_1 = avar_const_1,
    avar_const_2 = avar_const_2
  )
}

# ------------------------------------------------------------------------------
# Function: Asymptotic variance of the observed estimator Var(φ) on a dataset 
# generated by the functions pi1_0_fun, pi2_0_fun. 
# ------------------------------------------------------------------------------

compute_oif_variance <- function(data, estimates,
                                 pi1_0_fun, pi2_0_fun, 
                                 pi1_fun = NULL, pi2_fun = NULL, 
                                 pi1 = NULL, pi2 = NULL) {
  
  mom <- compute_moments(data, estimates, pi1_0_fun, pi2_0_fun)
  if (!is.null(pi1_fun)) {pi1 <- pi1_fun(cbind(data$X_B, data$X_t)) }
  if (!is.null(pi2_fun)) {pi2 <- pi2_fun(cbind(data$X_B, data$X_t, data$X_M))}
  
  
  # First expectation term: average over all samples
  term1_mean <- mean(mom$avar_const_1 * (mom$eps_t^2) / pi1)
  # Second expectation term: average over C >=2
  term2 <- (mom$eps_M %*% mom$avar_const_2)^2 / (pi1 * pi2)
  term2_mean <- weighted.mean(term2[data$C >= 2], w = mom$Weight_2inf[data$C >= 2])
  
  return(term1_mean + term2_mean)
}

# ------------------------------------------------------------------------------
# Function: Construct the control functions pi as a function of lambda
# ------------------------------------------------------------------------------

construct_pi_1 <- function(lambda, g1, g2_cond_exp, c1, c2_fun, 
                           indicator_2inf, weight_2inf, eps_M, n_sub,
                           xBt, x_M_base) {
  if (g1 >= lambda * c1) {
    return(1)
  } else {
    # select subset of n_sub eps_M samples to calculate the conditional expectation.
    n_sub <- 100
    n_2inf <- sum(indicator_2inf)
    n_sub <- min(n_sub, n_2inf)
    
    indicator_sub <- sort(sample(1:n_2inf, n_sub))
    indicator_full <- which(indicator_2inf)[indicator_sub]
    
    ones <- as.matrix(rep(1, n_sub))
    X_M_simulated <- (ones %*% x_M_base) + eps_M[indicator_full, ]
    X_B_simulated <- as.matrix(rep(1, n_sub)) %*% as.matrix(xBt$X_B)
    X_t_simulated <- rep(xBt$X_t, n_sub)
    data_simulated <-   data.frame(X_t = X_t_simulated)
    data_simulated$X_B <- I(X_B_simulated)
    data_simulated$X_M <- I(X_M_simulated)
    c2 <- c2_fun(data_simulated[, c("X_B", "X_t", "X_M")])
    
    #c2 <- c2_fun(cbind(as.matrix(rep(1, n_sub)) %*% as.matrix(xBt), X_M_simulated))
    
    weight_sub <- weight_2inf[indicator_full] * n_2inf / n_sub
    c2_cond_exp <- weighted.mean(c2, w = weight_sub)
    return(min(1, max(sqrt(g1 / (lambda * c1)), 
                      sqrt((g1 + g2_cond_exp) / (lambda * (c1 + c2_cond_exp)))
    )))
  }
}

construct_pi_2 <- function(lambda, g1, g2, c1, c2){
  return(ifelse(g1 >= lambda * c1,
                pmin(1, sqrt(g2 / (lambda * c2))),
                pmin(1, sqrt(g2 * c1 / (lambda * g1 * c2)))))
}

# ------------------------------------------------------------------------------
# Function: Compute the expected cost
# ------------------------------------------------------------------------------

# indicator: potential trimming of input vectors
compute_expected_cost <- function(pi1, pi2, c0, c1, c2, 
                                  weight_2inf, indicator = NULL) { 
  if (!is.null(indicator)) { 
    term <- c0 + (pi1[indicator] * c1[indicator]) + 
      (pi1[indicator] * pi2[indicator] * c2[indicator])
    weighted.mean(term, w = weight_2inf[indicator])
  } else {
    term <- c0 + (pi1 * c1) + (pi1 * pi2 * c2)
    weighted.mean(term, w = weight_2inf)
  }
}

# ------------------------------------------------------------------------------
# Function: Determine the optimal \pi = (\pi_1, \pi_2) under budget
# ------------------------------------------------------------------------------

compute_optimal_pi <- function(data, estimates,
                               c0 = 0,
                               c1_fun,
                               c2_fun,
                               b0,
                               pi1_0_fun, pi2_0_fun,
                               lambda_lower = 1e-8, lambda_upper = 1e8,
                               n_sub = 20000) {
  # Moments
  indicator_2inf <- data$C >= 2
  mom <- compute_moments(data, estimates,
                         pi1_0_fun, pi2_0_fun)
  X_M_base <- data$X_M - mom$eps_M
  g1 <- (mom$eps_t)^2 * mom$avar_const_1
  g2 <- (mom$eps_M %*% mom$avar_const_2)^2
  g2_cond_exp <- weighted.mean(g2, w = mom$Weight_2inf)
  c1 <- c1_fun(data[, c("X_B", "X_t")])
  c2 <- c2_fun(data[, c("X_B", "X_t", "X_M")])
  
  # compute expected cost
  compute_expected_cost_lambda <- function(lambda) {
    enumeration_2inf <- which(indicator_2inf)
    
    pi1 <- vapply(
      1:sum(indicator_2inf),
      function(i) {
        index_i <- enumeration_2inf[i]
        construct_pi_1(lambda, g1[index_i], g2_cond_exp, c1[index_i], c2_fun,
                       indicator_2inf, mom$Weight_2inf, mom$eps_M, n_sub,
                       data[index_i, c("X_B", "X_t")], X_M_base[index_i, ])}
      , numeric(1))
    pi2 <- construct_pi_2(lambda, g1, g2, c1, c2)[indicator_2inf]
    
    compute_expected_cost(pi1, pi2, c0, 
                          c1_fun(data[enumeration_2inf, c("X_B", "X_t")]),
                          c2_fun(data[enumeration_2inf, c("X_B", "X_t", "X_M")]),
                          mom$Weight_2inf[indicator_2inf])
  }
  
  # Solve for lambda
  fn_root <- function(lambda) compute_expected_cost_lambda(lambda) - b0
  lower <- lambda_lower
  upper <- lambda_upper
  root <- tryCatch({
    uniroot(fn_root, lower = lower, upper = upper)
  }, error = function(e) {
    # value to adapt for computational speed
    grid <- 10^seq(-2, 2, length.out = 200)
    vals <- sapply(grid, fn_root)
    idx <- which(diff(sign(vals)) != 0)
    if (length(idx) == 0) stop("Failed to bracket lambda. Adjust b0 or costs.")
    uniroot(fn_root, lower = grid[idx[1]], upper = grid[idx[1] + 1])
  })
  lambda_star <- root$root
  
  # evaluate optimal control on the given dataset
  pi1_star_fun <- function(data, lambda) {
    c1 <- c1_fun(data[, c("X_B", "X_t")])
    eps_t <- data$X_t - data$X_B %*% estimates$beta_tB
    g1 <- (eps_t)^2 * mom$avar_const_1
    
    X_M_base <- data$X_B %*% t(estimates$beta_MB) + data$X_t %*% estimates$beta_Mt
    eps_M_old <- mom$eps_M
    g2_old <- (eps_M_old %*% mom$avar_const_2)^2
    g2_cond_exp <- weighted.mean(g2_old, w = mom$Weight_2inf)
    
    pi1_star <- vapply(
      1:dim(data)[1],
      function(i) {
        construct_pi_1(lambda_star, g1[i], g2_cond_exp, c1[i], c2_fun,
                       indicator_2inf, mom$Weight_2inf, mom$eps_M, n_sub = sum(indicator_2inf),
                       data[i, c("X_B", "X_t")], X_M_base[i, ])}
      , numeric(1))
  }
  pi2_star_fun <- function(data, lambda_star) { 
    c1 <- c1_fun(data[, c("X_B", "X_t")])
    eps_t <- data$X_t - data$X_B %*% estimates$beta_tB
    g1 <- (eps_t)^2 * mom$avar_const_1
    eps_M_new <- data$X_M - data$X_t %*% estimates$beta_Mt - data$X_B %*% t(estimates$beta_MB)
    g2 <- (eps_M_new %*% mom$avar_const_2)^2
    c2 <- c2_fun(data[, c("X_B", "X_t", "X_M")])
    pi2_star <- construct_pi_2(lambda_star, g1, g2, c1, c2) 
  }
  
  list(
    lambda_star = lambda_star,
    pi12_star = cbind(pi1_star_fun(data, lambda_star), pi2_star_fun(data, lambda_star)),
    pi1_star_fun = pi1_star_fun,
    pi2_star_fun = pi2_star_fun
  )
}
