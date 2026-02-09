# ------------------------------------------------------------------------------
# This script tests the usefullness under misspecification
# ------------------------------------------------------------------------------

source("./Documents/papers/2025_sample_optimization/code/R/functions.R")

# ------------------------------------------------------------------------------
# Quadratic causal effect estimation
# ------------------------------------------------------------------------------

causal_effect_frontdoor <- function(data, pi1_fun, pi2_fun){
estimate_frontdoor <- function(data, t_grid, pi1_fun, pi2_fun) {
  pi1 <- pi1_fun(data[, c("X_B", "X_t")])# idx_MtB
  pi2 <- pi2_fun(data[, c("X_B", "X_t", "X_M")]) # idx_MtB
  
  data <- as.matrix(data)
  data <- as.data.frame(data)
  X_B <- data[, c("X_B.1", "X_B.2")]        
  X_M <- data[, c("X_M.1","X_M.2","X_M.3")]       
  X_t <- data[, "X_t", drop = FALSE]            
  X_r <- data[, "X_r", drop = FALSE]          
  
  ## Indices for available data
  idx_tB  <- which(!is.na(X_t))                 
  idx_MtB <- which(!is.na(X_M[,1]))            
  idx_rMtB <- which(!is.na(X_r))           
  
  
  ## 1) Model T | B
  df_tB <- data[idx_tB, c("X_B.1", "X_B.2", "X_t")]
  mod_t_given_B <- lm(X_t ~ ., data = df_tB)
  
  ## 2) Model M | T, B  (multivariate regression)
  df_MtB <- data[idx_MtB, c("X_B.1", "X_B.2", "X_t", "X_M.1","X_M.2","X_M.3")]
  df_MtB$X_t <- as.numeric(df_MtB$X_t)
  mod_M_given_tB <- lm(cbind(X_M.1, X_M.2, X_M.3) ~ X_B.1 + X_B.2 + poly(X_t, 3, raw = TRUE), 
                       data = df_MtB, weights = 1 / pi1[idx_MtB])
  
  ## 3) Model R | M, T, B
  df_rMtB <- data[idx_rMtB, c("X_B.1", "X_B.2", "X_t", "X_M.1","X_M.2","X_M.3", "X_r")]
  df_rMtB$X_t <- as.numeric(df_rMtB$X_t) 
  mod_r_given_MtB <- lm(X_r ~ ., data = df_rMtB, 
                        weights = 1 / (pi1[idx_rMtB] * pi2[idx_rMtB]))
  
  ## Baseline covariates for outer expectation
  B_mat <- X_B[idx_tB, , drop = FALSE]
  n_B   <- nrow(B_mat)
  
  ## Empirical distribution of T | B approximated by observed X_t
  T_vec <- X_t[idx_tB, , drop = FALSE]
  n_T   <- dim(T_vec)[1]
  
  ## Function computing E[R | do(T = t0)]
  E_r_do_t <- function(t0) {
    ER_given_B <- numeric(n_B)
    
    df_tB_copy <- df_tB
    df_tB_copy[, "X_t"] <- t0
    df_tB_copy$X_t <- as.numeric(df_tB_copy$X_t)
    m_hat_all <- predict(mod_M_given_tB, newdata = df_tB_copy)
    
    for (i in seq_len(n_B)) { #sample(seq_len(n_B), min(n_B, 2000))) {

      m_i <- m_hat_all[i, ]   # 1 × d_M
      
      ## Step 2: inner expectation over empirical T_j
      new_r <- as.data.frame(cbind(
        X_B = as.matrix(rep(1, n_T)) %*% as.matrix(B_mat[i, ]),
        X_t = as.matrix(X_t[idx_tB, , drop = FALSE]),
        X_M = matrix(rep(m_i, each = n_T), nrow = n_T)
      ))
      colnames(new_r)[4:6] <- c("X_M.1","X_M.2","X_M.3")
      
      r_hat_ij <- predict(mod_r_given_MtB, newdata = new_r)
      ER_given_B[i] <- mean(r_hat_ij)
    }
    
    mean(ER_given_B[ER_given_B!=0])
  }
  
  ## Compute E[R | do(t)] on grid
  E_vals <- sapply(t_grid, E_r_do_t)
  
  list(t_grid = t_grid,
       E_r_do_t = E_vals)
}

t_grid <- seq(-0.1, 0.1, length.out = 10)
fd_res <- estimate_frontdoor(data, t_grid, pi1_fun, pi2_fun) 
fd_res$E_r_do_t
}


# ------------------------------------------------------------------------------
# Generate Quadratic Data
# ------------------------------------------------------------------------------

generate_data_quadratic <- function(n, dB = 2, dM = 3, dS = 1,
                          beta_tB, beta_MB, beta_Mt,
                          beta_rB, beta_rM, beta_Mt2,
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun, pi2_fun) {
  
  X_B <- eps_B_fun(n)
  eps_tr <- eps_tr_fun(n)
  X_t <- X_B %*% beta_tB + eps_tr[, 1]
  eps_M <- eps_M_fun(n)
  X_M <-  X_B %*% t(beta_MB) + X_t %*% beta_Mt + X_t^2 %*% beta_Mt2 + eps_M
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
# Experiment
# ------------------------------------------------------------------------------

set.seed(123)

dB = 2
dM = 3
dS = 1
beta_tB = c(0.5, -0.2)
beta_MB = rbind(c(0.3, 0.1), c(0.5, 0.2), c(-0.1, 0.3))
beta_Mt = c(0.7, 0.2, 0.1)
beta_Mt2 = c(0.1, 0.2, 0.4)
beta_rB = c(0.2, -0.1)
beta_rM = c(0.5, 0.4, -0.3)

sigma_B = rbind(c(1, 0.7), c(0.7, 1.5))
eps_B_fun <- function(n) {mvtnorm::rmvt(n, sigma = sigma_B, df = 5)}
sigma_tr <- rbind(c(1, -0.5), c(-0.5, 1.5))
eps_tr_fun <- function(n) {eps_tr0 <- mvtnorm::rmvt(n, sigma = sigma_tr, df = 5)}
sigma_M = rbind(c(1, 0.3, 0), c(0.3, 1.5, -0.5), c(0, -0.5, 1))
eps_M_fun <- function(n) {MASS::mvrnorm(n, rep(0, 3), sigma_M)}

c1_fun <- function(XtB) 0.1 * apply(XtB, 1, function(xtB) sqrt(sum(xtB^2))) 
c2_fun <- function(XBtM) { rep(0.5, nrow(XBtM))} 

c0 <- 0

# ------------------------------------------------------------------------------

# generate data and optimize
pi1_fun <- function(XBt) { rep(1, nrow(XBt)) }
pi2_fun <- function(XBtM) { ifelse(apply(XBtM, MARGIN = 1, function(x) {any(is.na(x))}), 
                                   NA, 1) }

set.seed(123)

data <- generate_data_quadratic(1000, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM,beta_Mt2,
                      eps_B_fun, eps_tr_fun, eps_M_fun,
                      pi1_fun,
                      pi2_fun)
estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
indicator_2inf <- data$C >= 2
mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)

# evaluate vanilla
avar_vanilla <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                     pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                      c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                      mom$Weight_2inf, indicator = indicator_2inf)
print(c("asymptotic variance:", avar_vanilla))
print(c("the associated average cost is:", cost_vanilla))
scale <- 1.5
b0 <- cost_vanilla / scale # scale for equal cost

# optimize for the control function
opt <- compute_optimal_pi(data, estimates, 
                          c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                          b0 = b0, 
                          pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
avar_opt <- compute_oif_variance(data, estimates,
                                 pi1_fun, pi2_fun,
                                 pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
cost_opt <- compute_expected_cost(opt$pi12_star[, 1], opt$pi12_star[, 2],
                                  c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                  mom$Weight_2inf, indicator = indicator_2inf)

n_reps <- 50
sample_sizes <- c(100, 250, 500, 750, 1000, 2500, 5000)#, 7500)#, 10000)

# sample and estimate on vanilla
set.seed(123)
xi_hat_vanilla <- array(NA, dim = c(n_reps, length(sample_sizes), 10))
for (i in 1:n_reps){
  for (j in 1:length(sample_sizes)) {
    print(Sys.time())
    print(c("vanilla", i, j))
    data <- generate_data_quadratic(sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, beta_Mt2,
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    xi_hat_vanilla[i, j,] <- causal_effect_frontdoor(data, pi1_fun, pi2_fun)
  }
}

# sample and estimate on optimized
set.seed(123)
pi1_star_fun <- function(data) opt$pi1_star_fun(data, opt$lambda_star)
pi2_star_fun <- function(data) opt$pi2_star_fun(data, opt$lambda_star)
xi_hat_opt <- array(NA, dim = c(n_reps, length(sample_sizes), 10))
for (i in 1:n_reps){
  for (j in 1:length(sample_sizes)) {
    print(Sys.time())
    print(c("opt", i, j))
    data <- generate_data_quadratic(scale * sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, beta_Mt2,
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_star_fun,
                          pi2_star_fun)
    xi_hat_opt[i, j, ] <- causal_effect_frontdoor(data, pi1_star_fun, pi2_star_fun)
  }
}

# mse_computation:
t_grid <- seq(-0.1, 0.1, length.out = 10) 
xi <- sum(beta_Mt * beta_rM) * t_grid + sum(beta_Mt2 * beta_rM) * t_grid^2 

mse_vanilla <- array(NA, dim = length(sample_sizes))
mse_opt <- array(NA, dim = length(sample_sizes))
for (s in 1:length(sample_sizes)) {
  mse_vanilla[s] <- ModelMetrics::mse(array(xi_hat_vanilla[, s, ]), 
                                         rep(xi, each = n_reps))
  mse_opt[s] <- ModelMetrics::mse(array(xi_hat_opt[, s, ]), 
                                         rep(xi, each = n_reps))
}
#
saveRDS(list(xi_hat_opt = xi_hat_opt, xi_hat_vanilla = xi_hat_vanilla,
             mse_opt = mse_opt, mse_vanilla = mse_vanilla, 
             sample_sizes = sample_sizes, 
             lambda_star = opt$lambda_star, 
             avar_vanilla = avar_vanilla, 
             avar_opt = avar_opt, 
             pi12_star = opt$pi12_star),
        "./Documents/papers/2025_sample_optimization/code/R/runs/calibration/main_miss.RDS")

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/calibration_miss.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0)) 
plot(sample_sizes[1:7], mse_vanilla[1:7], log = "xy", col = "blue", type = "l", 
     main = "Model Misspecification: Non-Linear Data",
     xlab = expression("n ·" ~ b[0]),
     ylab = expression(MSE(hat(xi)[n])),
     ylim = c(min(c(mse_opt, mse_vanilla)), max(c(mse_opt, mse_vanilla))),
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
mtext(paste("Avg Cost Full-Data =", scale, "· Avg Cost Observed-Data"), side = 3, line = 0.5, cex = 0.8)
lines(sample_sizes[1:7], mse_opt[1:7], col = "red", lwd = 3)
legend("topright", 
       legend = c("Full Data", "Optimized", "Simulated"), 
       col = c("blue", "red", "black"), 
       pch    = c(19, 19, NA),   
       lty    = c(NA, NA, 1),
       lwd    = c(NA, NA, 2.5),
       cex = 1.2)
grid()
dev.off()
