source("./Documents/papers/2025_sample_optimization/code/R/function_back_door.R")

# ------------------------------------------------------------------------------
# Set up 

set.seed(123)
n = 1000
dB = 2
dM = 3
dS = 1

beta_tB = c(0.5, -0.2)
beta_MB = rbind(c(0.3, 0.1), c(0.5, 0.2), c(-0.1, 0.3))
beta_Mt = c(0.7, 0.2, 0.1)

sigma_B = rbind(c(1, 0.7), c(0.7, 1.5))
eps_B_fun <- function(n) {mvtnorm::rmvt(n, sigma = sigma_B, df = 5)}
sigma_tr <- rbind(c(1, -0.5), c(-0.5, 1.5))
eps_tr_fun <- function(n) {mvtnorm::rmvt(n, sigma = sigma_tr, df = 5)}
sigma_M = rbind(c(1, 0.3, 0), c(0.3, 1.5, -0.5), c(0, -0.5, 1))
eps_M_fun <- function(n) {MASS::mvrnorm(n, rep(0, 3), sigma_M)}

pi1_fun <- function(XBt) { rep(0.3, nrow(XBt)) + runif(nrow(XBt), -0.1, 0.1) }
pi2_fun <- function(XBtM) { ifelse(apply(XBtM, MARGIN = 1, 
                                         function(x) {any(is.na(x))}), NA, 
                                   rep(0.5, nrow(XBtM)) + 
                                     runif(nrow(XBtM), -0.1, 0.1)) }

# ------------------------------------------------------------------------------
# Generate data and estimate parameters

data <- generate_data(n, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                      eps_B_fun, eps_tr_fun, eps_M_fun,
                      pi1_fun)
estimates <- estimate_parameters(data, pi1_fun)
print(estimates)

# ------------------------------------------------------------------------------
# Find an optimal control function and compare performances

# Define costs:
c1_fun <- function(XtB) 0.1 * apply(XtB, 1, function(xtB) sqrt(sum(xtB^2))) #rep(0.1, nrow(XBt))
c0 <- 0

# base analysis
indicator_2inf <- data$C >= 2
mom <- compute_moments(data, estimates, pi1_fun)
avar_base <- compute_oif_variance(data, estimates,
                                  pi1_fun,
                                  pi1_fun = pi1_fun)
cost_base <- compute_expected_cost(pi1_fun(data),
                                   c0, c1_fun(data[, c("X_B", "X_t")]), 
                                   mom$Weight_2inf, indicator = indicator_2inf)
print(c("asymptotic variance:", avar_base))
print(c("the associated average cost is:", cost_base))
b0 <- cost_base

cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]),
                                      c0, c1_fun(data[, c("X_B", "X_t")]), 
                                      mom$Weight_2inf, indicator = indicator_2inf)
avar_vanilla <- compute_oif_variance(data, estimates, pi1_fun,
                                     pi1 = rep(1, dim(data)[1]))

# optimize the control function
start.time <- Sys.time()
opt <- compute_optimal_pi(data, estimates, 
                          c0 = 0.0, c1_fun = c1_fun,
                          b0 = b0, 
                          pi1_0_fun = pi1_fun, n_sub = 10000)
Sys.time() - start.time

print(c("optimal lambda is", opt$lambda))
print("optimal pi head is:"); print(head(opt$pi12_star))

avar_opt <- compute_oif_variance(data, estimates,
                                 pi1_fun, 
                                 pi1 = opt$pi12_star[, 1])
cost_opt <- compute_expected_cost(opt$pi1_star, 
                                  c0, c1_fun(data[, c("X_B", "X_t")]),
                                  mom$Weight_2inf, indicator = indicator_2inf)
print(c("asymptotic variance:", avar_opt))
print(c("the associated average cost is:", cost_opt))

