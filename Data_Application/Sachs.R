library(readr)
# load data
Sachs <- read_csv("Documents/papers/2025_sample_optimization/code/R/data_analysis/Sachs.csv")
Sachs <- as.data.frame(Sachs)
n <- nrow(Sachs)
# load functions
source("./Documents/papers/2025_sample_optimization/code/R/functions.R")

pi1_fun <- function(XBt) { rep(1, nrow(XBt)) }
pi2_fun <- function(XBtM) { ifelse(apply(XBtM, MARGIN = 1, function(x) {any(is.na(x))}), 
                                   NA, 1) }
indicator_2inf <- rep(TRUE, n)


# Create different front door scenarios
# B = PKA, PKC, t = Raf, M = Mek, r = Erk
# B = PKA, PKC, t = Raf, M = Mek, r = Akt
# B = PKA, PKC, t = Raf, M = Erk, r = Akt: scale = 1.3500000 -> rel eff = 0.7661884
# B = PKA,      t = Mek, M = Erk, r = Akt: scale = 1.5400000 -> rel eff = 0.6801851

names <- colnames(Sachs)
Sachs <- as.matrix(Sachs)
colnames(Sachs) <- NULL

Sachs1 <- data.frame(C = rep(Inf, n), X_t = Sachs[, 1], X_M = Sachs[, 2], X_r = Sachs[, 6])
X_B <- Sachs[, c(8, 9)]
Sachs1$X_B = I(X_B)

Sachs2 <- data.frame(C = rep(Inf, n), X_t = Sachs[, 1], X_M = Sachs[, 2], X_r = Sachs[, 7])
X_B <- Sachs[, c(8, 9)]
Sachs2$X_B = I(X_B)


Sachs3 <- data.frame(C = rep(Inf, n), X_t = Sachs[, 1], X_M = Sachs[, 6], X_r = Sachs[, 7])
X_B <- Sachs[, c(8, 9)]
Sachs3$X_B = I(X_B)

Sachs4 <- data.frame(C = rep(Inf, n), X_t = Sachs[, 2], X_M = Sachs[, 6], X_r = Sachs[, 7])
X_B <- Sachs[, c(8)]
Sachs4$X_B = I(X_B)

# ------------------------------------------------------------------------------
# Sachs1
# ------------------------------------------------------------------------------

c1_fun <-  function(XBt) { rep(1, nrow(XBt))}
c2_fun <- function(XBtM) { rep(1, nrow(XBtM))} 

# ------------------------------------------------------------------------------
# Sachs1
# ------------------------------------------------------------------------------

set.seed(1234)
data <- Sachs1
estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)

# evaluate vanilla
avars_vanilla <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                      pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                      c0 = 3, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                      mom$Weight_2inf, indicator = indicator_2inf)
for (scale in seq(1.001, 1.01, by = 0.001)) {
  b0 <- cost_vanilla / scale
  opt <- compute_optimal_pi(data, estimates, 
                            c0 = 3, c1_fun = c1_fun, c2_fun = c2_fun,
                            b0 = b0, 
                            pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
  avars_opt <- compute_oif_variance(data, estimates,
                                    pi1_fun, pi2_fun,
                                    pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  print(c(scale, avars_opt / avars_vanilla / scale))
}

# no: smaller 1.001000

# ------------------------------------------------------------------------------
# Sachs2
# ------------------------------------------------------------------------------

set.seed(1234)
data <- Sachs2
estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)

# evaluate vanilla
avars_vanilla <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                      pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                      c0 = 3, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                      mom$Weight_2inf, indicator = indicator_2inf)
for (scale in seq(1.001, 2, by = 0.1)) {
  b0 <- cost_vanilla / scale # scale for equal cost
  opt <- compute_optimal_pi(data, estimates, 
                            c0 = 3, c1_fun = c1_fun, c2_fun = c2_fun,
                            b0 = b0, 
                            pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
  avars_opt <- compute_oif_variance(data, estimates,
                                    pi1_fun, pi2_fun,
                                    pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  print(c(scale, avars_opt / avars_vanilla / scale))
}
# no: smaller 1.0010

# ------------------------------------------------------------------------------
# Sachs3 -> no optimum found
# ------------------------------------------------------------------------------

set.seed(1234)
data <- Sachs3
estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)

# evaluate vanilla
avars_vanilla <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                      pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                      c0 = 3, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                      mom$Weight_2inf, indicator = indicator_2inf)
for (scale in seq(1.35, 1.4, by = 0.01)) {
  b0 <- cost_vanilla / scale # scale for equal cost
  opt <- compute_optimal_pi(data, estimates, 
                            c0 = 3, c1_fun = c1_fun, c2_fun = c2_fun,
                            b0 = b0, 
                            pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
  avars_opt <- compute_oif_variance(data, estimates,
                                    pi1_fun, pi2_fun,
                                    pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  print(c(scale, avars_opt / avars_vanilla / scale))
}
# scale = [1] 1.3500000 0.7661884

# ------------------------------------------------------------------------------
# Sachs4 
# ------------------------------------------------------------------------------

set.seed(1234)
data <- Sachs4
estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)

# evaluate vanilla
avars_vanilla <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                      pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                      c0 = 2, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                      mom$Weight_2inf, indicator = indicator_2inf)
for (scale in seq(1.4, 1.6, by = 0.01)) {
  b0 <- cost_vanilla / scale # scale for equal cost
  opt <- compute_optimal_pi(data, estimates, 
                            c0 = 2, c1_fun = c1_fun, c2_fun = c2_fun,
                            b0 = b0, 
                            pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
  avars_opt <- compute_oif_variance(data, estimates,
                                    pi1_fun, pi2_fun,
                                    pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  print(c(scale, avars_opt / avars_vanilla / scale))
}
# scale = [1] 1.5400000 0.6801851