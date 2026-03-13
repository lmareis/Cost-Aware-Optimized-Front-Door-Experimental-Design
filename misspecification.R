# ------------------------------------------------------------------------------
# This script tests the usefullness under misspecification
# ------------------------------------------------------------------------------

source("./Documents/papers/2025_sample_optimization/code/R/functions.R")

# ------------------------------------------------------------------------------
# Quadratic causal effect estimation
# ------------------------------------------------------------------------------

estimate_M <- function(data, beta_tB, pi1_fun, quadratic = FALSE) {
  eps_t <- data$X_t - as.matrix(data$X_B) %*% beta_tB
  W <- ifelse(data$C >= 2, 1 / pi1_fun(data[, c("X_B", "X_t")]), 0)
  
  if (quadratic) {
    fit <- lm(as.matrix(data$X_M) ~ data$X_t + I(data$X_t^2) + as.matrix(data$X_B),
              weights = W)
    if (!is.matrix(coef(fit))) {
      return(list(
        beta_Mt  = as.matrix(coef(fit)[2]),
        beta_Mt2 = as.matrix(coef(fit)[3]),
        beta_MB  = t(coef(fit)[-c(1, 2, 3)])
      ))
    }
    return(list(
      beta_Mt  = coef(fit)[2, ],
      beta_Mt2 = coef(fit)[3, ],
      beta_MB  = t(coef(fit)[-c(1, 2, 3), ])
    ))
  } else {
    fit <- lm(as.matrix(data$X_M) ~ data$X_t + as.matrix(data$X_B),
              weights = W)
    if (!is.matrix(coef(fit))) {
      return(list(
        beta_Mt  = as.matrix(coef(fit)[2]),
        beta_Mt2 = NULL,
        beta_MB  = t(coef(fit)[-c(1, 2)])
      ))
    }
    return(list(
      beta_Mt  = coef(fit)[2, ],
      beta_Mt2 = NULL,
      beta_MB  = t(coef(fit)[-c(1, 2), ])
    ))
  }
}

estimate_parameters <- function(data, pi1_fun, pi2_fun, quadratic = FALSE) {
  beta_tB <- estimate_tB(data)
  M_est   <- estimate_M(data, beta_tB, pi1_fun, quadratic = quadratic)
  r_est   <- estimate_r(data, beta_tB, M_est$beta_Mt, M_est$beta_MB,
                        pi1_fun, pi2_fun)
  
  list(beta_tB  = beta_tB,
       beta_Mt  = M_est$beta_Mt,
       beta_Mt2 = M_est$beta_Mt2,
       beta_MB  = M_est$beta_MB,
       beta_rB  = as.matrix(r_est[1:length(beta_tB)]),
       gamma    = r_est[length(beta_tB) + 1],
       beta_rM  = as.matrix(r_est[-(1:(length(beta_tB) + 1))]))
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
set.seed(127)

dB = 2
dM = 3
dS = 1
beta_tB = c(0.5, -0.2)
beta_MB = rbind(c(0.3, 0.1), c(0.5, 0.2), c(-0.1, 0.3))
beta_Mt = c(0.7, 0.2, 0.1)
beta_Mt2 = -c(0.1, 0.2, 0.4)
beta_rB = c(0.2, -0.1)
beta_rM = c(0.5, 0.4, -0.3)
t_grid <- seq(-0.1, 0.1, length.out = 10) 


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


data <- generate_data_quadratic(1000, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM,beta_Mt2,
                      eps_B_fun, eps_tr_fun, eps_M_fun,
                      pi1_fun,
                      pi2_fun)
estimates <- estimate_parameters(data, pi1_fun, pi2_fun, quadratic = TRUE)
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
sample_sizes <- c(100, 250, 500, 750, 1000, 2500, 5000, 7500)
n_sizes <- length(sample_sizes)

# sample and estimate on vanilla
xi_hat_vanilla <- array(NA, dim = c(n_reps, length(sample_sizes), 10))
cost_vanilla <- array(NA, dim = c(n_reps, length(sample_sizes)))
for (i in 1:n_reps){
  for (j in 1:length(sample_sizes)) {
    print(Sys.time())
    print(c("vanilla", i, j))
    data <- generate_data_quadratic(sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, beta_Mt2,
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun, quadratic = TRUE)
    xi_hat_vanilla[i, j, ] <- sum(estimates$beta_Mt * estimates$beta_rM) * t_grid +
      sum(estimates$beta_Mt2 * estimates$beta_rM) * t_grid^2
    cost_vanilla[i, j] <- sum(c1_fun(data[, c("X_B", "X_t")]) * (data[, "C"] > 1) + c2_fun(data) * (data[, "C"] > 2))
  }
}

# sample and estimate on optimized
set.seed(123)
pi1_star_fun <- function(data) opt$pi1_star_fun(data, opt$lambda_star)
pi2_star_fun <- function(data) opt$pi2_star_fun(data, opt$lambda_star)
xi_hat_opt <- array(NA, dim = c(n_reps, length(sample_sizes), 10))
cost_opt <- array(NA, dim = c(n_reps, length(sample_sizes)))
for (i in 1:n_reps){
  for (j in 1:length(sample_sizes)) {
    print(Sys.time())
    print(c("opt", i, j))
    data <- generate_data_quadratic(scale * sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, beta_Mt2,
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_star_fun,
                          pi2_star_fun)
    cost_opt[i, j] <- sum(c1_fun(data[, c("X_B", "X_t")]) * (data[, "C"] > 1) + c2_fun(data) * (data[, "C"] > 2))
    estimates <- estimate_parameters(data, pi1_star_fun, pi2_star_fun, quadratic = TRUE)
    xi_hat_opt[i, j, ] <- sum(estimates$beta_Mt * estimates$beta_rM) * t_grid +
      sum(estimates$beta_Mt2 * estimates$beta_rM) * t_grid^2
  }
}

# mse_computation:
xi <- sum(beta_Mt * beta_rM) * t_grid + sum(beta_Mt2 * beta_rM) * t_grid^2 
mse_vanilla <- array(NA, dim = length(sample_sizes))
mse_opt <- array(NA, dim = length(sample_sizes))
for (s in 1:length(sample_sizes)) {
  mse_vanilla[s] <- ModelMetrics::mse(array(xi_hat_vanilla[, s, ]), 
                                         rep(xi, each = n_reps))
  mse_opt[s] <- ModelMetrics::mse(array(xi_hat_opt[, s, ]), 
                                         rep(xi, each = n_reps))
}

saveRDS(list(xi_hat_opt = xi_hat_opt, xi_hat_vanilla = xi_hat_vanilla,
             mse_opt = mse_opt, mse_vanilla = mse_vanilla, 
             sample_sizes = sample_sizes, 
             lambda_star = opt$lambda_star, 
             avar_vanilla = avar_vanilla, 
             avar_opt = avar_opt, 
             pi12_star = opt$pi12_star),
        "./Documents/papers/2025_sample_optimization/code/R/runs/calibration/main_miss.RDS")

RDS <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/calibration/main_miss.RDS")
xi_hat_opt <- RDS$xi_hat_opt
xi_hat_vanilla <- RDS$xi_hat_vanilla
mse_opt <- RDS$mse_opt
mse_vanilla <- RDS$mse_vanilla
sample_sizes <- RDS$sample_sizes
lambda_star <- RDS$lambda_star
avar_vanilla <- RDS$avar_vanilla
avar_opt <- RDS$avar_opt
pi12_star <- RDS$pi12_star

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
lines(sample_sizes[1:7] * colMeans(cost_opt)[1:7] / colMeans(cost_vanilla)[1:7], mse_opt[1:7], col = "red", lwd = 3)
legend("topright", 
       legend = c("Full Data", "Optimized", "Simulated"), 
       col = c("blue", "red", "black"), 
       pch    = c(19, 19, NA),   
       lty    = c(NA, NA, 1),
       lwd    = c(NA, NA, 2.5),
       cex = 1.2)
grid()
dev.off()
