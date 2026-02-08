# ------------------------------------------------------------------------------
# This script evaluates the time and samples to fix the optimal lambda
# ------------------------------------------------------------------------------

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

pi1_fun <- function(XBt) { rep(1, nrow(XBt)) }
pi2_fun <- function(XBtM) { ifelse(apply(XBtM, MARGIN = 1, function(x) {any(is.na(x))}), 
                                   NA, 1) }

set.seed(123)

data <- generate_data(1000, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
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
b0 <- cost_vanilla / scale # half the cost

# optimize for the control function

set.seed(123)

n_reps <- 50
sample_sizes <- c(100, 250, 500, 750, 1000, 2500, 5000, 10000)
compute_times <- array(NA, dim = c(n_reps, length(sample_sizes)))
avars_opt <- array(NA, dim = c(n_reps, length(sample_sizes)))
lambda_stars <- array(NA, dim = c(n_reps, length(sample_sizes)))
costs_opt <- array(NA, dim = c(n_reps, length(sample_sizes)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(sample_sizes)) {
    print(j)
    data <- generate_data(sample_sizes[j], dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    start.time <- Sys.time()
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    compute_times[i,j] <- as.numeric(Sys.time() - start.time, units = "secs")
    lambda_stars[i,j] <- opt$lambda_star
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                     pi1_fun, pi2_fun,
                                     pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
    costs_opt[i,j] <- compute_expected_cost(opt$pi12_star[, 1], opt$pi12_star[, 2],
                                      c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data[, c("X_B", "X_t", "X_M")]),
                                      mom$Weight_2inf, indicator = data$C >= 2)
  }
}

saveRDS(list(var_opt = var_opt,
             sample_sizes = sample_sizes, 
             compute_times = compute_times, 
             lambda_stars = lambda_stars,
             avars_opt = avars_opt,
             costs_opt = costs_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/computational_sensitivity.RDS")
RDS_lambda <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/computational_sensitivity.RDS")
var_opt <- RDS_lambda$var_opt
sample_sizes <- RDS_lambda$sample_sizes
compute_times <- RDS_lambda$compute_times
lambda_stars <- RDS_lambda$lambda_stars
avars_opt <- RDS_lambda$avars_opt
costs_opt <- RDS_lambda$costs_opt

# Compute means and variances across replicates
means <- list(
  compute_times = colMeans(compute_times, na.rm = TRUE),
  avars_opt     = colMeans(avars_opt, na.rm = TRUE),
  lambda_stars  = colMeans(lambda_stars, na.rm = TRUE),
  costs_opt     = colMeans(costs_opt, na.rm = TRUE)
)

vars <- list(
  compute_times = apply(compute_times, 2, var, na.rm = TRUE),
  avars_opt     = apply(avars_opt, 2, var, na.rm = TRUE),
  lambda_stars  = apply(lambda_stars, 2, var, na.rm = TRUE),
  costs_opt     = apply(costs_opt, 2, var, na.rm = TRUE)
)

# Plot means with error bars showing variance
par(mfrow = c(2,2))  # 2x2 grid of plots

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/computational_sensitivity/ct.pdf",
    width = 8, height = 6)
par(mgp = c(3, 1, 0), mar = c(5, 5, 4, 2))
plot(sample_sizes, means$compute_times, type="b", log = "xy", pch=19, col="blue",
     ylim = range(c(means$compute_times - sqrt(vars$compute_times), 
                    means$compute_times, 
                    means$compute_times + sqrt(vars$compute_times))),
     main="Compute Time (s)", xlab="Sample Size", ylab="Mean ± SD",
     cex.main=3, cex.lab=2.5, cex.axis = 2, lwd = 3)
arrows(sample_sizes, means$compute_times - sqrt(vars$compute_times),
       sample_sizes, means$compute_times + sqrt(vars$compute_times),
       angle=90, code=3, length=0.05, col="blue")
dev.off()

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/computational_sensitivity/avar.pdf",
    width = 8, height = 6)
par(mgp = c(3, 1, 0), mar = c(5, 5, 4, 2))
plot(sample_sizes, means$avars_opt, type="b", log = "xy", pch=19, col="red",
     ylim = range(c(means$avars_opt - sqrt(vars$avars_opt), 
                    means$avars_opt, 
                    means$avars_opt + sqrt(vars$avars_opt))),
     main="Computed Asymp Variance", xlab="Sample Size", ylab="Mean ± SD",
     cex.main=3, cex.lab=2.5, cex.axis = 2, lwd = 3)
arrows(sample_sizes, means$avars_opt - sqrt(vars$avars_opt),
       sample_sizes, means$avars_opt + sqrt(vars$avars_opt),
       angle=90, code=3, length=0.05, col="red")
dev.off()

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/computational_sensitivity/lambda_star.pdf",
    width = 8, height = 6)
par(mgp = c(3, 1, 0), mar = c(5, 5, 4, 2))
plot(sample_sizes, means$lambda_stars, type="b", log = "xy", pch=19, col="darkgreen",
     ylim = range(c(means$lambda_stars - sqrt(vars$lambda_stars), 
                    means$lambda_stars, 
                    means$lambda_stars + sqrt(vars$lambda_stars))),
     main=expression(bold("Computed Lambda Star" ~ lambda^"*")), xlab="Sample Size", ylab="Mean ± SD",
     cex.main=3, cex.lab=2.5, cex.axis = 2, lwd = 3)
arrows(sample_sizes, means$lambda_stars - sqrt(vars$lambda_stars),
       sample_sizes, means$lambda_stars + sqrt(vars$lambda_stars),
       angle=90, code=3, length=0.05, col="darkgreen")
dev.off()

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/computational_sensitivity/avg_cost.pdf",
    width = 8, height = 6)
par(mgp = c(3, 1, 0), mar = c(5, 5, 4, 2))
plot(sample_sizes, means$costs_opt, type="b", log = "xy", pch=19, col="purple",
     ylim = range(c(means$costs_opt - sqrt(vars$costs_opt), 
                    means$costs_opt, 
                    means$costs_opt + sqrt(vars$costs_opt))),
     main="Computed Avg Cost", xlab="Sample Size", ylab="Mean ± SD",
     cex.main=3, cex.lab=2.5, cex.axis = 2, lwd = 3)
arrows(sample_sizes, means$costs_opt - sqrt(vars$costs_opt),
       sample_sizes, means$costs_opt + sqrt(vars$costs_opt),
       angle=90, code=3, length=0.05, col="purple")
dev.off()




