# File Purpose: This R script is designed to parallelizing a single run of a monte carlo simulation
# In the case of debugging and ensuring that everything is working specifically

library(tidyverse)
library(MASS)
library(boot)
library(extraDistr)
library(lme4)
library(parallel)

########################################

# Parameters to change!
num_iter <- 500
num_boot <- 100
num_clusters <- 30
cluster_range <- c(100, 300)
ICC <- 0.01
independent <- F # Is classification independent of covariates (T = Yes, F = No)
parametric <- T # Should probabilities be calculated parametrically? (T = Yes/Model 1, F = No/Model 2)
randomized <- F # Are we in a randomized or an observational context?
clustered <- T # Dependent on if we are in an observational or randomized setting

#######################################

calculate_ATE_prop <- function(data, parametric, clustered, randomized){
  validation <- data |>
    filter(!is.na(Y_ij_0) | !is.na(Y_ij_1))
  
  # Fitting a propensity score
  ## Model assumes perfectly IID observations
  if(randomized){
      pi_hat <- mean(data$T_ij)
  } else{
    if(!clustered){
      pi_hat_model <- glm(T_ij ~ x1 + x2 + x3 + x4, data = data, family = "binomial")
      pi_hat <- predict(pi_hat_model, newdata = data, type = 'response')
    } else{
      pi_hat_model <- glmer(T_ij ~ x1 + x2 + x3 + x4 + (1 | cluster_num), data = data, family = "binomial")
      pi_hat <- predict(pi_hat_model, newdata = data, type = 'response')
    }
  }
  
  # Fitting Probability models, of the form P(Y_ij, T_ij, X_ij)
  ## Removing previously computed probabilities
  data <- data |>
    dplyr::select(-starts_with("p"))
  
  if (parametric) {
    expit <- glm(
      Y_ij_star ~ (T_ij * Y_ij) + (T_ij * x1) + (T_ij * x2) + (T_ij *x3) + (x4),
      data = validation,
      family = "binomial"
    )
    
    ## Fitting P(0, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 0, T_ij = 1)
    p_0_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 1, T_ij = 1)
    p_1_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(0, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 0, T_ij = 0)
    p_0_0 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 1, T_ij = 0)
    p_1_0 <- predict(expit, newdata = temp, type = 'response')
  } else{
    p_0_1 <- sum(validation$Y_ij_star * (1 - validation$Y_ij) * validation$T_ij) / sum((1 - validation$Y_ij) * validation$T_ij)
    
    p_1_1 <- sum(validation$Y_ij_star * validation$Y_ij * validation$T_ij) / sum(validation$Y_ij * validation$T_ij)
    
    p_0_0 <- sum(validation$Y_ij_star * (1 - validation$Y_ij) * (1 - validation$T_ij)) / sum((1 - validation$Y_ij) * (1 - validation$T_ij))
    
    p_1_0 <- sum(validation$Y_ij_star * validation$Y_ij * (1 - validation$T_ij)) / sum(validation$Y_ij * (1 - validation$T_ij))
  }
  
  realistic_data <- cbind(data, p_0_1, p_1_1, p_0_0, p_1_0)
  
  # Calculating the actual SSW ATE
  
  E_1 <- sum(((realistic_data$T_ij * realistic_data$Y_ij_star) - (pi_hat * realistic_data$p_0_1)
  ) / (pi_hat * (
    realistic_data$p_1_1 - realistic_data$p_0_1
  ))) / nrow(realistic_data)
  E_0 <- sum((((1 - realistic_data$T_ij) * realistic_data$Y_ij_star
  ) - ((1 - pi_hat) * realistic_data$p_0_0)) / ((1 - pi_hat) * (
    realistic_data$p_1_0 - realistic_data$p_0_0))) / nrow(realistic_data)
  return(E_1 - E_0)
}

boot_ATE <- function(data, parametric, clustered, randomized){
  sample_ATE <- replicate(num_boot, {
    # Should this be by individual or by clusters? Seemingly Clusters
    sample_idx <- sample(1:length(unique(data$cluster_num)), size = 30 ,replace = T)
    boot_data <- lapply(seq_along(sample_idx), function(idx){
      data[data$cluster_num == sample_idx[idx], ] |>
        mutate(cluster_num = idx)
    })
    boot_data <- data.table::rbindlist(boot_data)
    calculate_ATE_prop(boot_data, parametric, clustered, randomized)
  })
  return(c(quantile(sample_ATE, probs = 0.025), quantile(sample_ATE, probs = 0.975)))
}

ATE <- rep(NA, num_iter)
ATE_est <- rep(NA, num_iter)
coverage <- rep(NA, num_iter)

one_loop <- function(i){
  data <- list()
  for (j in 1:num_clusters) {
    # Generating our covariate matrix X_i per cluster
    cluster_size <- rdunif(1, cluster_range[1], cluster_range[2])
    x1 <- rnorm(cluster_size, mean = 1, sd = 1)
    x2 <- mvrnorm(mu = rep(0.5, cluster_size), Sigma = ((0.5 * diag(x = cluster_size)) + (
      0.05 * matrix(
        data = 1,
        nrow = cluster_size,
        ncol = cluster_size
      )
    )))
    x3 <- rbinom(cluster_size, 1, 0.55)
    x4 <- rep(runif(1), cluster_size)
    cluster_level_data <- rbind(x1, x2, x3, x4) #X_i
    # Generating our observed outcome, Y_ij, Y*_ij, and V_ij
    # Starting with Y_ij(0) and Y_ij(1)
    b_i_y <- rnorm(1, 0, sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
    Y_ij_0_pi <- inv.logit((-1) + (c(0.15, 0.2, 0.15, -0.15) %*% cluster_level_data) + b_i_y) # Bernoulli parameter for calculating Y_ij(0)
    Y_ij_0 <- rbern(length(Y_ij_0_pi), Y_ij_0_pi)
    Y_ij_1_pi <- inv.logit((-0.25) + (c(0.15, 0.2, 0.15, -0.15) %*% cluster_level_data) + b_i_y) # Bernoulli parameter for calculating Y_ij(1)
    Y_ij_1 <- rbern(length(Y_ij_1_pi), Y_ij_1_pi)
    # Now moving onto Y*_ij(0) and Y*_ij(1)
    b_i_ystar <- 0 # READ: Paper stated that they kept this at zero?
    if(independent){
      Y_ij_0_star_pi <-  inv.logit(-1.25 + (1.5 * Y_ij_0) + b_i_ystar)# Bernoulli parameter for calculating Y*_ij(0)
      Y_ij_star_0 <- rbern(length(Y_ij_0_star_pi), Y_ij_0_star_pi)
      Y_ij_1_star_pi <-  inv.logit(-1 + (2.5 * Y_ij_1) + b_i_ystar) # Bernoulli parameter for calculating Y*_ij(1)
      Y_ij_star_1 <- rbern(length(Y_ij_1_star_pi), Y_ij_1_star_pi)
    } else{
      Y_ij_0_star_pi <-  inv.logit(-1.75 + (c(0.25, -0.25, -0.15, -0.1) %*% cluster_level_data) + (1.5 * Y_ij_0) + b_i_ystar)# Bernoulli parameter for calculating Y*_ij(0)
      Y_ij_star_0 <- rbern(length(Y_ij_0_star_pi), Y_ij_0_star_pi)
      Y_ij_1_star_pi <-  inv.logit(-1.25 + (c(-0.25, -0.15, -0.25, -0.1) %*% cluster_level_data) + (2.5 * Y_ij_1) + b_i_ystar) # Bernoulli parameter for calculating Y*_ij(1)
      Y_ij_star_1 <- rbern(length(Y_ij_1_star_pi), Y_ij_1_star_pi)
    }
    # Now, we calculate V_ij
    b_iv <- rnorm(1, 0, sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
    V_ij_0_pi <- inv.logit((-0.25) + (c(-0.5, -0.5, 0.25, -0.25) %*% cluster_level_data) + (-0.15 * Y_ij_0)  + b_iv) # Bernoulli parameter for calculating V_ij(0)
    V_ij_0 <- rbern(length(V_ij_0_pi), V_ij_0_pi)
    V_ij_1_pi <- inv.logit((-0.5) + (c(-0.5, -0.5, 0.25, -0.25) %*% cluster_level_data) + (0.15 * Y_ij_1)  + b_iv)  # Bernoulli parameter for calculating V_ij(1)
    V_ij_1 <- rbern(length(V_ij_1_pi), V_ij_1_pi)
    # Finally, We randomly assign our treatment (READ: this is the biggest difference between the first file and this one
    # assume IID for now, even if this logically doesn't make sense in a cluster based setting)
    if(randomized){
      T_ij <- rbern(1)
    } else{
      T_ij_random <- rnorm(1, mean = 1, sd = 2)
      T_ij_pi <- inv.logit(0.5 + (c(0.15, 0.25, -0.1, 0.1) %*% cluster_level_data) + T_ij_random)
      T_ij <- rbern(length(T_ij_pi), T_ij_pi)
    }
    cluster_num <- j
    cluster_level_data <- rbind(
      cluster_level_data,
      Y_ij_0,
      Y_ij_1,
      Y_ij_star_0,
      Y_ij_star_1,
      V_ij_0,
      V_ij_1,
      T_ij,
      cluster_num
    )
    data[[j]] <- t(cluster_level_data)
  }
  true_data <- do.call(rbind, lapply(data, as.data.frame))
  # Calculating the true value of the ATE
  ATE <- mean(true_data$Y_ij_1) - mean(true_data$Y_ij_0)
  # Estimating the SSW ATE
  ## Creating the realistic data set, accounting for the potential outcomes and missingness of gold standard
  realistic_data <- true_data |>
    mutate(
      V_ij_0 = if_else(T_ij == 0, V_ij_0, NA),
      V_ij_1 = if_else(T_ij == 1, V_ij_1, NA),
      Y_ij_star_0 = if_else(T_ij == 0, Y_ij_star_0, NA),
      Y_ij_star_1 = if_else(T_ij == 1, Y_ij_star_1, NA),
      Y_ij_star = coalesce(Y_ij_star_0, Y_ij_star_1),
      Y_ij_0 = if_else(V_ij_0 == 1, Y_ij_0, NA),
      Y_ij_1 = if_else(V_ij_1 == 1, Y_ij_1, NA),
      Y_ij = coalesce(Y_ij_0, Y_ij_1)
    )
  
  ATE_est <- calculate_ATE_prop(realistic_data, parametric, clustered, randomized)
  CI_boot <- boot_ATE(realistic_data, parametric, clustered, randomized)
  # CI_boot returns c(Lower, Upper) bootstrapped confidence intervals
  coverage <- CI_boot[1] < ATE & ATE < CI_boot[2]
  return(list(ATE = ATE, ATE_est = ATE_est, coverage = coverage))
}


RNGkind("L'Ecuyer-CMRG")
set.seed(999)
system.time(final_data <- mclapply(1:num_iter, one_loop, mc.cores = 50))
final_data <- do.call(rbind, lapply(final_data, as.data.frame))
write_csv(final_data, "Single-Results.csv")
