# ------------------------------------------------------------------------------
# This script compares the optimized against the full estimation and variance
# ------------------------------------------------------------------------------

# Setup
source("./Documents/papers/2025_sample_optimization/code/R/functions.R")

set.seed(123)

dB = 2
dM = 3
dS = 1
beta_tB = c(0.5, -0.2)
beta_MB = rbind(c(0.3, 0.1), c(0.5, 0.2), c(-0.1, 0.3))
beta_Mt = c(0.7, 0.2, 0.1)
beta_rB = c(0.2, -0.1)
beta_rM = c(0.5, 0.4, -0.3)
xi <- sum(beta_Mt * beta_rM)

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

data <- generate_data(1000, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM,  # TODO: perhaps redo with more base samples
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
sample_sizes <- c(100, 250, 500, 750, 1000, 2500, 5000, 7500)#, 10000)

# sample and estimate on vanilla
set.seed(123)
xi_hat_vanilla <- array(NA, dim = c(n_reps, length(sample_sizes)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(sample_sizes)) {
    data <- generate_data(sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    xi_hat_vanilla[i, j] <- sum(estimates$beta_Mt * estimates$beta_rM) - xi
  }
}
var_vanilla <- apply(xi_hat_vanilla+ xi, 2, function(xi_hat) ModelMetrics::mse(xi_hat, rep(xi, n_reps)))

# sample and estimate on optimized
set.seed(123)
pi1_star_fun <- function(data) opt$pi1_star_fun(data, opt$lambda_star)
pi2_star_fun <- function(data) opt$pi2_star_fun(data, opt$lambda_star)
xi_hat_opt <- array(NA, dim = c(n_reps, length(sample_sizes)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(sample_sizes)) {
    print(j)
    data <- generate_data(scale * sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_star_fun,
                          pi2_star_fun)
    # probably best to generate on full data and coarsen afterwards.
    estimates <- estimate_parameters(data, pi1_star_fun, pi2_star_fun)
    xi_hat_opt[i, j] <- sum(estimates$beta_Mt * estimates$beta_rM) - xi
  }
}
var_opt <- apply(xi_hat_opt+ xi, 2, function(xi_hat) ModelMetrics::mse(xi_hat, rep(xi, n_reps)))
saveRDS(list(var_opt = var_opt, var_vanilla = var_vanilla, 
             sample_sizes = sample_sizes, 
             lambda_star = opt$lambda_star, 
             avar_vanilla = avar_vanilla, 
             avar_opt = avar_opt, 
             pi12_star = opt$pi12_star),
        "./Documents/papers/2025_sample_optimization/code/R/runs/calibration/main.RDS")

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/calibration.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0)) 
plot(sample_sizes[1:7], var_vanilla[1:7], log = "xy", col = "blue", type = "l", 
     main = "Empirical and Theoretical MSE Across Budget Levels",
     xlab = expression("n ·" ~ b[0]),
     ylab = expression(MSE(hat(xi)[n]) ~~ "and" ~~ MSE[asymp](hat(xi)[n])),
     ylim = c(min(c(var_opt, var_vanilla)), max(c(var_opt, var_vanilla))),
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
mtext(paste("Avg Cost Full-Data =", scale, "· Avg Cost Observed-Data"), side = 3, line = 0.5, cex = 0.8)
lines(sample_sizes[1:7], var_opt[1:7], col = "red", lwd = 3)
legend("topright", 
       legend = c("Full Data", "Optimized", "Simulated", "Theory"), 
       col = c("blue", "red", "black", "black"), 
       pch    = c(19, 19, NA, NA),   
       lty    = c(NA, NA, 1, 2),
       lwd    = c(NA, NA, 2.5, 2.5),
       cex = 1.2)

# full data sqrt{n} (\hat \xi_n - \xi) \sim N(0, avar_vanilla)
# -> n Var(\hat \xi_n) = avar_vanilla
# Var(\hat \xi_n) = avar_vanilla / n
var_vanilla_theory <- avar_vanilla / sample_sizes
lines(sample_sizes[1:7], var_vanilla_theory[1:7], col = "blue", lty = 2.5)
var_opt_theory <- avar_opt / sample_sizes / scale
lines(sample_sizes[1:7], var_opt_theory[1:7], col = "red", lty = 2.5)
grid()
dev.off()

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/calibration_pi_star_distribution.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0)) 
plot(opt$pi12_star[, 1], opt$pi12_star[, 2], pch = 1,
     main = expression(bold("Optimal Propensity" ~ pi^"*" ~ "Distribution")),
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5,
     xlab = expression(pi[1]^"*"), 
     ylab = expression(pi[2]^"*"), 
)
dev.off()



  