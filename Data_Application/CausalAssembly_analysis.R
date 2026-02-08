library(readr)
library(igraph)
station1 <- read_csv("Documents/papers/2025_sample_optimization/code/R/data_analysis/station1.csv")
source("./Documents/papers/2025_sample_optimization/code/R/function_back_door.R")

# ------------------------------------------------------------------------------
# Estimate DAG

data_mat <- as.matrix(station1)
suffStat <- list(C = cor(data_mat), n = nrow(data_mat))

pc_fit <- pc(
  suffStat = suffStat,
  indepTest = gaussCItest,   
  alpha = 0.03,             
  labels = colnames(station1),
  verbose = TRUE
)
summary(pc_fit)
plot(pc_fit, main = "Estimated DAG with PC Algorithm")


# ------------------------------------------------------------------------------
# CausalAssembly
# ------------------------------------------------------------------------------

n <- nrow(station1)
c1_fun <-  function(XBt) { rep(1, nrow(XBt))}
pi1_fun <- function(XBt) { rep(1, nrow(XBt)) }
indicator_2inf <- rep(TRUE, n)

station1 <- as.matrix(station1)
colnames(station1) <- NULL
CA <- data.frame(C = rep(Inf, n), X_t = station1[, 6], X_M = station1[, 5], X_B = station1[, 4])

set.seed(1234)
data <- CA
estimates <- estimate_parameters(data, pi1_fun)
mom <- compute_moments(data, estimates, pi1_fun)

# evaluate vanilla
avars_vanilla <- compute_oif_variance(data, estimates, pi1_fun, 
                                      pi1 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]),
                                      c0 = 2, c1_fun(data[, c("X_B", "X_t")]),
                                      mom$Weight_2inf, indicator = indicator_2inf)
for (scale in seq(1.01, 1.1, by = 0.01)){
  print(scale)
  b0 <- cost_vanilla / scale
  opt <- compute_optimal_pi(data, estimates, 
                            c0 = 2, c1_fun = c1_fun,
                            b0 = b0, 
                            pi1_0_fun = pi1_fun, n_sub = 10000)
  avars_opt <- compute_oif_variance(data, estimates,
                                    pi1_fun,
                                    pi1 = opt$pi1_star)
  print(c(scale, avars_opt / avars_vanilla / scale))
  
}

# 1.0900000 -> 0.9469889
