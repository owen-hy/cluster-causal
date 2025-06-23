# File Purpose: This R script is designed to account for an observational setting
# Specifically, we are now dealing with a scenario where {Y_ij(0), Y_ij(1)} _|_ T_ij | X_ij
# Here we also assume that treatment is assigned at the cluster level

library(tidyverse)
library(MASS)
library(boot)
library(extraDistr)
library(nnet)
library(parallel)

num_iter <- 30
num_clusters <- 30

calculate_ATE_prop <- function(data, parametric, clustered){
  validation <- data |>
    filter(!is.na(Y_ij_0) | !is.na(Y_ij_1))
  # Fitting a propensity score
  # OY: TODO
  propensity_ind <- multinom(cluster_num ~ x1 + x2 + x3, data = data)
  pi_hat 
  
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

boot_ATE <- function(data, parametric, clustered){
  sample_ATE <- replicate(100, {
    # Should this be by individual or by clusters? Seemingly Clusters
    sample_idx <- sample(1:length(unique(data$cluster_num)), size = 30 ,replace = T)
    boot_data <- lapply(seq_along(sample_idx), function(idx){
      data[data$cluster_num == sample_idx[idx], ] |>
        mutate(cluster_num = idx)
    })
    boot_data <- data.table::rbindlist(boot_data)
    calculate_ATE_prop(boot_data, parametric, clustered)
  })
  return(c(quantile(sample_ATE, probs = 0.025), quantile(sample_ATE, probs = 0.975)))
}

ATE_sim_one_obs <- function(cluster_range, parametric, ICC, independent, clustered){
  ATE <- rep(NA, num_iter)
  ATE_est <- rep(NA, num_iter)
  coverage <- rep(NA, num_iter)
  for(i in 1:num_iter) {
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
      # Finally, We randomly assign our treatment 
      # OY: TODO
      x5 <- rbinom(1, 1, 0.5)# Adding one more cluster level covariate
      T_ij_pi <- inv.logit((-0.1 * x4) + (0.25 * x5))
      T_ij <- rbinom(1, 1, T_ij_pi)
      cluster_num <- j
      cluster_level_data <- rbind(
        cluster_level_data,
        x5,
        Y_ij_0,
        Y_ij_1,
        Y_ij_star_0,
        Y_ij_star_1,
        V_ij_0,
        V_ij_1,
        T_ij,
        as.factor(cluster_num)
      )
      data[[j]] <- t(cluster_level_data)
    }
    true_data <- do.call(rbind, lapply(data, as.data.frame))
    # Calculating the true value of the ATE
    ATE[i] <- mean(true_data$Y_ij_1) - mean(true_data$Y_ij_0)
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
    
    ATE_est[i] <- calculate_ATE_prop(realistic_data, parametric, clustered)
    CI_boot <- boot_ATE(realistic_data, parametric, clustered)
    # CI_boot returns c(Lower, Upper) bootstrapped confidence intervals
    coverage[i] <- CI_boot[1] < ATE[i] & ATE[i] < CI_boot[2]
  }
  return(list(ATE = mean(ATE), ATE_est = mean(ATE_est), 
              variance = var(ATE_est),
              coverage = mean(coverage),
              bias_se = sqrt(var(ATE_est) / num_iter),
              cov_se = sqrt((mean(coverage) * (1 - mean(coverage))) / num_iter),
              size = if_else(cluster_range[1] == 100, "small", "large"),
              parametric = parametric, ICC = ICC, independent = independent))
}


