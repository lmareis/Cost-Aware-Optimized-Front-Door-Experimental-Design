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
# beta_tB # ~ 25 min execution time
# ------------------------------------------------------------------------------

set.seed(122)

n_reps <- 50
scales <- c(0.01, 0.05, 0.1, 0.15, 1/ 5:1, 2:5, 10, 20)
start.time <- Sys.time()
avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    data <- generate_data(500, dB, dM, dS, scales[j] * beta_tB, beta_MB, beta_Mt, beta_rB, beta_rM, 
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
             avars_opt = avars_opt,
             costs_opt = costs_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_tB.RDS")
RDS_tB <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_tB.RDS")
avars_opt <- RDS_tB$avars_opt
avars_rel <- RDS_tB$avars_rel
scales <- RDS_tB$scales

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/beta_tB.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression(beta[tB] ~ "Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Parameter" ~ beta[tB])),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="blue",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="blue")

# add right axis
axis(side=4, col="blue", col.axis="blue",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="blue", cex=1.5)

dev.off()
Sys.time() - start.time


# ------------------------------------------------------------------------------
# beta_MB
# ------------------------------------------------------------------------------

set.seed(123)

n_reps <- 50
scales <- c(0.01, 0.05, 0.1, 0.15, 1/ 5:1, 2:5, 10, 20)

avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    data <- generate_data(500, dB, dM, dS, beta_tB, scales[j] * beta_MB, beta_Mt, beta_rB, beta_rM, 
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
             avars_opt = avars_opt,
             costs_opt = costs_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_MB.RDS")
RDS_MB <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_MB.RDS")
avars_opt <- RDS_MB$avars_opt
avars_rel <- RDS_MB$avars_rel
scales <- RDS_MB$scales

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/beta_MB.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression(beta[MB] ~ "Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Parameter" ~ beta[MB])),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="blue",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="blue")

# add right axis
axis(side=4, col="blue", col.axis="blue",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="blue", cex=1.5)

dev.off()


# ------------------------------------------------------------------------------
# beta_Mt
# ------------------------------------------------------------------------------

set.seed(123)

n_reps <- 50
scales <- c(0.01, 0.05, 0.1, 0.15, 1/ 5:1, 2:5, 10, 20)

avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, scales[j] * beta_Mt, beta_rB, beta_rM, 
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
             avars_opt = avars_opt,
             costs_opt = costs_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_Mt.RDS")
RDS_Mt <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_Mt.RDS")
avars_opt <- RDS_Mt$avars_opt
avars_rel <- RDS_Mt$avars_rel
scales <- RDS_Mt$scales

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/beta_Mt.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression(beta[Mt] ~ "Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Parameter" ~ beta[Mt])),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="blue",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="blue")

# add right axis
axis(side=4, col="blue", col.axis="blue",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="blue", cex=1.5)

dev.off()


# ------------------------------------------------------------------------------
# beta_rB
# ------------------------------------------------------------------------------

set.seed(124)

n_reps <- 50
scales <- c(0.01, 0.05, 0.1, 0.15, 1/ 5:1, 2:5, 10, 20)

avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, scales[j] * beta_rB, beta_rM, 
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
             avars_opt = avars_opt,
             costs_opt = costs_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_rB.RDS")
RDS_rB <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_rB.RDS")
avars_opt <- RDS_rB$avars_opt
avars_rel <- RDS_rB$avars_rel
scales <- RDS_rB$scales

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/beta_rB.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression(beta[rB] ~ "Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Parameter" ~ beta[rB])),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="blue",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="blue")

# add right axis
axis(side=4, col="blue", col.axis="blue",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="blue", cex=1.5)

dev.off()


# ------------------------------------------------------------------------------
# beta_rM
# ------------------------------------------------------------------------------

set.seed(123)

n_reps <- 50
scales <- c(0.01, 0.05, 0.1, 0.15, 1/ 5:1, 2:5, 10, 20)

avars_vanilla <- array(NA, dim = c(n_reps, length(scales)))
avars_opt <- array(NA, dim = c(n_reps, length(scales)))
for (i in 1:n_reps){
  print(i)
  for (j in 1:length(scales)) {
    print(j)
    data <- generate_data(500, dB, dM, dS, beta_tB, beta_MB, beta_Mt, beta_rB, scales[j] * beta_rM, 
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
             avars_opt = avars_opt,
             costs_opt = costs_opt),
        "./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_rM.RDS")
RDS_rM <- readRDS("./Documents/papers/2025_sample_optimization/code/R/runs/parameter_sensitivity_rM.RDS")
avars_opt <- RDS_rM$avars_opt
avars_rel <- RDS_rM$avars_rel
scales <- RDS_rM$scales

var_opt <- apply(avars_opt, 2, var)
mean_opt <- apply(avars_opt, 2, mean)
var_rel <- apply(avars_rel, 2, var)
mean_rel <- apply(avars_rel, 2, mean)
sd_opt <- sqrt(var_opt)
sd_rel <- sqrt(var_rel)

pdf("./Documents/papers/2025_sample_optimization/code/R/plots/parameter_sensitivity/beta_rM.pdf",
    width = 8, height = 6)
par(mgp = c(2, 0.5, 0), mar = c(5, 5, 4, 6)) 

# left axis: opt, log-log
plot(scales, mean_opt, type="b", pch=19, col="black",
     ylim=range(mean_opt - sd_opt, mean_opt + sd_opt),
     xlab=expression(beta[rM] ~ "Scaling"), ylab="Optimized Variance",
     log="xy", main = expression(bold("Parameter" ~ beta[rM])),
     cex.main=2, cex.lab=1.5, cex.axis = 1.5, lwd = 3)
grid()

# error bars for opt
arrows(scales, mean_opt - sd_opt, scales, mean_opt + sd_opt,
       angle=90, code=3, length=0.05, col="black")

# overlay rel on right axis, log-log
par(new=TRUE)
plot(scales, mean_rel, type="b", pch=19, col="blue",
     ylim=range(mean_rel - sd_rel, mean_rel + sd_rel),
     axes=FALSE, xlab="", ylab="", log="xy",
     cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, lwd = 3)

# error bars for rel
arrows(scales, mean_rel - sd_rel, scales, mean_rel + sd_rel,
       angle=90, code=3, length=0.05, col="blue")

# add right axis
axis(side=4, col="blue", col.axis="blue",
     cex.lab=1.5, cex.axis = 1.5)
mtext("Relaltive Efficiency", side=4, line=3, col="blue", cex=1.5)

dev.off()


