# ------------------------------------------------------------------------------
# This script evaluates the sensitivity of optimization w.r.t model pars
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

scale <- 1.5 # cost scale

# ------------------------------------------------------------------------------

# generate data and optimize
pi1_fun <- function(XBt) { rep(1, nrow(XBt)) }
pi2_fun <- function(XBtM) { ifelse(apply(XBtM, MARGIN = 1, function(x) {any(is.na(x))}), 
                                   NA, 1) }

# ------------------------------------------------------------------------------
# sigma_B # ~ 15 min execution time
# ------------------------------------------------------------------------------

set.seed(122)

n_reps <- 50
scales <- seq(0, sqrt(3), length.out = 11)
start.time <- Sys.time()
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
cors <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_B = rbind(c(1, scales[j] * 0.7), c(scales[j] * 0.7, 1.5))
    eps_B_fun <- function(n) {mvtnorm::rmvt(n, sigma = sigma_B, df = 5)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
    cors[i, j] <- cor(eps_B_fun(1000))[1,2]
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine


saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt,
             cors = cors),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_B.RDS")
RDS_B <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_B.RDS")
avars_rel <- RDS_B$avars_rel
cors <- RDS_B$cors
avars_opt <- RDS_B$avars_opt

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
corB <- apply(cors, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_B.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(corB, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Cor(" ~ epsilon[B1]~ "," ~ epsilon[B2]~ ")"), ylab="Optimized Variance",
     log="y", main = expression(bold("Nuisance Cor(" ~ epsilon[B1] ~ "," ~ epsilon[B2] ~ ")")),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(corB, mean_opt - sd_opt, corB, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(corB, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="y",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(corB, mean_rel - sd_rel, corB, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()
Sys.time() - start.time


# ------------------------------------------------------------------------------
# sigma_varB # ~ 15 min execution time
# ------------------------------------------------------------------------------

set.seed(123)

n_reps <- 50
scales <- c(0.6, 0.75, 1, 1.3, 1.75, 2.25, 3, 4.5, 7, 10, 15)
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
cors <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_B = rbind(c(scales[j] * 1, 0.7), c(0.7, scales[j] * 1.5))
    eps_B_fun <- function(n) {mvtnorm::rmvt(n, sigma = sigma_B, df = 5)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine

saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_varB.RDS")
RDS_varB <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_varB.RDS")
avars_rel <- RDS_varB$avars_rel
scales <- RDS_varB$scales
avars_opt <- RDS_varB$avars_opt

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_varB.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Var(" ~ epsilon[B1]~ "), and Var(" ~ epsilon[B2]~ ") Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Nuisance Var(" ~ epsilon[B1] ~ "), and Var(" ~ epsilon[B2] ~ ")")),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()

# ------------------------------------------------------------------------------
# sigma_t #
# ------------------------------------------------------------------------------

set.seed(124)

n_reps <- 50
scales <- c(0.2, 0.4, 0.7, 1, 1.5, 2.5, 5, 8, 12.5, 20)
start.time <- Sys.time()
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_tr <- rbind(c(scales[j] * 1, -0.5), c(-0.5, 1.5))
    eps_tr_fun <- function(n) {eps_tr0 <- mvtnorm::rmvt(n, sigma = sigma_tr, df = 5)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine

saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_t.RDS")
RDS_t <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_t.RDS")
avars_rel <- RDS_t$avars_rel
scales <- RDS_t$scales
avars_opt <- RDS_t$avars_opt

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_t.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Var(" ~ epsilon[t]~ ") Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Nuisance Var(" ~ epsilon[t] ~ ")" )),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()

# ------------------------------------------------------------------------------
# sigma_r #
# ------------------------------------------------------------------------------

set.seed(125)

n_reps <- 50
scales <- c(0.3, 0.5, 0.7, 1, 1.5, 2.5, 5, 8, 12.5, 20)
start.time <- Sys.time()
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_tr <- rbind(c(1, -0.5), c(-0.5, scales[j] *1.5))
    eps_tr_fun <- function(n) {eps_tr0 <- mvtnorm::rmvt(n, sigma = sigma_tr, df = 5)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine

saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_r.RDS")
RDS_r <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_r.RDS")
avars_rel <- RDS_r$avars_rel
scales <- RDS_r$scales
avars_opt <- RDS_r$avars_opt

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)

sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_r.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Var(" ~ epsilon[r]~ ") Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Nuisance Var(" ~ epsilon[r] ~ ")" )),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()


# ------------------------------------------------------------------------------
# sigma_tr #
# ------------------------------------------------------------------------------

set.seed(126)

n_reps <- 50
scales <- -seq(0, 1.2, length.out = 11)
start.time <- Sys.time()
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
cors <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_tr <- rbind(c(1, scales[j] * (-0.5)), c(scales[j] * (-0.5), 1.5))
    eps_tr_fun <- function(n) {eps_tr0 <- mvtnorm::rmvt(n, sigma = sigma_tr, df = 5)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
    cors[i, j] <- cor(eps_tr_fun(1000))[1,2]
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine

saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt,
             cors = cors),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_tr.RDS")
RDS_tr <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_tr.RDS")
avars_rel <- RDS_tr$avars_rel
avars_opt <- RDS_tr$avars_opt
cors <- RDS_tr$cors

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
cortr <- apply(cors, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_tr.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(cortr, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Cor(" ~ epsilon[t] ~ "," ~ epsilon[r] ~ ")"), ylab="Optimized Variance",
     log="y", main = expression(bold("Nuisance Cor(" ~ epsilon[t] ~ "," ~ epsilon[r] ~ ")" )),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(cortr, mean_opt - sd_opt, cortr, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(cortr, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="y",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(cortr, mean_rel - sd_rel, cortr, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()

# ------------------------------------------------------------------------------
# sigma_varM # ~ 15 min execution time
# ------------------------------------------------------------------------------

set.seed(127)

n_reps <- 50
scales <- c(0.5, 0.6, 0.8, 1.1, 1.6, 2.25, 3, 4.5, 7, 10, 15)
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
cors <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_M = rbind(c(scales[j] *1, 0.3, 0), c(0.3, scales[j] *1.5, -0.5), c(0, -0.5, scales[j] *1))
    eps_M_fun <- function(n) {MASS::mvrnorm(n, rep(0, 3), sigma_M)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine


saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_varM.RDS")
RDS_varM <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_varM.RDS")
avars_rel <- RDS_varM$avars_rel
scales <- RDS_varM$scales
avars_opt <- RDS_varM$avars_opt

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_varM.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Var(" ~ epsilon[M1]~ "), and Var(" ~ epsilon[M2]~ "), and Var(" ~ epsilon[M3] ~ ") Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Nuisance Var(" ~ epsilon[M1] ~ "), Var(" ~ epsilon[M2] ~ "), and Var(" ~ epsilon[M3] ~ ")")),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()

# ------------------------------------------------------------------------------
# sigma_corM # ~ 15 min execution time
# ------------------------------------------------------------------------------

set.seed(128)

n_reps <- 50
scales <- seq(0, 2, by = 0.2)
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
cors <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    sigma_M = rbind(c(1, scales[j] *0.3, 0), c(scales[j] *0.3, 1.5, scales[j] *(-0.5)), c(0, scales[j] *(-0.5), 1))
    eps_M_fun <- function(n) {MASS::mvrnorm(n, rep(0, 3), sigma_M)}
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
                          eps_B_fun, eps_tr_fun, eps_M_fun,
                          pi1_fun,
                          pi2_fun)
    estimates <- estimate_parameters(data, pi1_fun, pi2_fun)
    indicator_2inf <- data$C >= 2
    mom <- compute_moments(data, estimates, pi1_fun, pi2_fun)
    
    # evaluate vanilla
    avars_vanilla[i,j] <- compute_oif_variance(data, estimates, pi1_fun, pi2_fun,
                                               pi1 = rep(1, dim(data)[1]), pi2 = rep(1, dim(data)[1]))
    cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]), rep(1, dim(data)[1]),
                                          c0, c1_fun(data[, c("X_B", "X_t")]), c2_fun(data),
                                          mom$Weight_2inf, indicator = indicator_2inf)
    b0 <- cost_vanilla / scale # scale for equal cost
    
    # optimize for the control function
    opt <- compute_optimal_pi(data, estimates, 
                              c0 = 0.0, c1_fun = c1_fun, c2_fun = c2_fun,
                              b0 = b0, 
                              pi1_0_fun = pi1_fun, pi2_0_fun = pi2_fun, n_sub = 10000)
    avars_opt[i,j] <- compute_oif_variance(data, estimates,
                                           pi1_fun, pi2_fun,
                                           pi1 = opt$pi12_star[, 1], pi2 = opt$pi12_star[, 2])
    cor_M <- cor(eps_M_fun(1000))
    cors[i,j] <- mean(c(abs(cor_M[1,2]), abs(cor_M[2,3])))
  }
}

avars_rel <- avars_opt / avars_vanilla / 1.5 # should be fine
var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
corM <- apply(cors, 2, mean)

saveRDS(list(avars_rel = avars_rel,
             scales = scales, 
             n_reps = n_reps,
             avars_opt = avars_opt,
             cors = cors),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_corM.RDS")
RDS_corM <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_corM.RDS")
avars_rel <- RDS_corM$avars_rel
cors <- RDS_corM$cors
avars_opt <- RDS_corM$avars_opt

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
corM <- apply(cors, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/sigma_corM.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(corM, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression("Cor(" ~ epsilon[M1]~ "," ~ epsilon[M2]~ "), and Cor(" ~ epsilon[M2]~ "," ~ epsilon[M3] ~ ")"), ylab="Optimized Variance",
     log="y", main = expression(bold("Nuisance Cor(" ~ epsilon[M1]~ "," ~ epsilon[M2]~ "), and Cor(" ~ epsilon[M2]~ "," ~ epsilon[M3] ~ ")")),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(corM, mean_opt - sd_opt, corM, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(corM, mean_rel, type="b", pch=19, col="red",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="y",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(corM, mean_rel - sd_rel, corM, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="red")

# add right axis
axis(side=4, col="red", col.axis="red",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="red", cex=1.5)

dev.off()

