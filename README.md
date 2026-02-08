# Cost-Aware-Optimized-Front-Door-Experimental-Design

This repository contains the code to the paper "Cost-Aware Optimized Front-Door Experimental Design" by (anonymized).

## Files:

|  |  |
|----|----|
| experiment_calibration.R | Code for Figure 2 / Supplement 3 |
| experiment_lambda_star.R | Code for Figure 3 |
| experiment_sensitivity.R | Code for blue plots in Figure 4 |
| experiment_sensitivity_nuisance.R | Code for red plots in Figure 4 |
| functions.R | Implementations of all front-door methods |
| misspecification.R | Code for Supplement 4 |

## Folders:

|  |  |
|------------------------------------|------------------------------------|
| Back_Door | Code for Back-door methods and experiments |
| Data_Application | Code for real-life applications. Datasets are not included due to licensing |
| plots | Plot of all Figures in the paper |
| runs | .RDS files of all runs used for the visualizations |

## Functions.R

| Method | Usage |
|------------------------------------|------------------------------------|
| generate_data | Generates data according to the observed-data multivariate additive noise linear front-door model $\mathcal{M}_\pi$ |
| estimate_tB | Estimates the linear coefficient $\beta_{tB}$ |
| estimate_M | Estimates the linear coefficients $\beta_{Mt}$ and $\beta_{MB}$. |
| estimate_r | Estimates the linear coefficients $\beta_{rM}$ and $\beta_{rB}$. |
| estimate_parameters | Estimates all coefficients in \$ \beta \$ as well as the best linear approximation \$ \gamma \$ of \$ E[\\varepsilon_r \| \varepsilon\_t] \$ |
| compute_moments | Computes residuals \$ \varepsilon\_t \$ and $\varepsilon_M$, the levarages \$ g_1 \$ and \$ g_2 \$ from Thm. 11 (called avar_1 and avar_2), and the propensity vector (Weight_2inf) |
| compute_oif_variance | Computes the observed-data efficient influence function |
| construct_pi_1 | Constructs \$\\pi_1\^\\ast \$ from Thm. 11 |
| construct_pi_2 | Constructs \$\\pi_2\^\\ast \$ from Thm. 11 |
| compute_expected_cost | Computes the expected cost of a sample under model \$ \mathcal{M}\_{\\pi} \$. |
| compute_optimal_pi | Computes the optimized \$ \pi\^\\ast \$ satisfying the budges constraint for a \$ \lambda\^\\ast \$. |

| Argument | Interpretation |
|------------------------------------|------------------------------------|
| n | Sample Size |
| dB | Dimension of \$ X_B \$ |
| dM | Dimension of \$ X_M \$ |
| dS | Dimension of \$ X_S \$ |
| beta_tB | Linear Coefficient \$ \beta\_{tB} \$ |
| beta_MB | Linear Coefficient \$ \beta\_{MB} \$ |
| beta_Mt | Linear Coefficient \$ \beta\_{Mt} \$ |
| beta_rB | Linear Coefficient \$ \beta\_{rB} \$ |
| beta_rM | Linear Coefficient \$ \beta\_{rM} \$ |
| eps_B_fun | Error function \$ \varepsilon\_B \$ |
| eps_tr_fun | Error function \$ \varepsilon\_{tr} \$ |
| eps_M_fun | Error function \$ \varepsilon\_M \$ |
| pi1_fun | Propensity function \$ \pi\_1 \$ |
| pi2_fun | Propensity function \$ \pi\_1 \$ |
| data | Dataset |
| estimates | A list list(beta_tB, beta_Mt, beta_MB, beta_rB, gamma, beta_rM) of estimates |
| pi1_0_fun | Initial propensity function \$ \pi\_1 \$ |
| pi2_0_fun | Initial propensity function \$ \pi\_2 \$ |
| pi1 | Vector of optimized propensities \$ \pi\_1 \$ |
| pi2 | Vector of optimized propensities \$ \pi\_2 \$ |
| c0 | Base cost of \$ X\_{Bt} \$ |
| c1_fun | Cost function \$ c_1(X\_{Bt}) \$ for \$ X_M \$ |
| c2_fun | Cost function \$ c_1(X\_{BtM}) \$ for \$ X_r \$ |
| b0 | Budget |
| lambda_lower | Lower bound for \$ \lambda \$ in the computational optimization |
| lambda_upper | Upper bound for \$ \lambda \$ in the computational optimization |
| n_sub | Sample size for Monte-Carlo sub-sampling of \$ E[g_2 \| \varepsilon\_{Bt} ] \$ |
