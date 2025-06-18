# File Purpose: This R script is designed to replicate the result of Dane's estimator in a simulation study
# This is done under the ideal situation where all clusters are randomized (Thus no propensity model is needed)

library(tidyverse)
library(MASS)
library(boot)
library(extraDistr)
library(geepack)
library(parallel)

num_iter <- 500 # Dane used 5000, just doing 100 until I learn how to parallel code
num_clusters <- 30

calculate_ATE <- function(data, parametric){
  validation <- data |>
    filter(!is.na(Y_ij_0) | !is.na(Y_ij_1))
  pi_hat <- mean(data$treatment)
  
  # Fitting Probability models, of the form P(Y_ij, T_ij, X_ij)
  
  ## Removing previously computed probabilities
  data <- data |>
    dplyr::select(-starts_with("p"))
  
  if (parametric) {
    # Fitting logistic model
    expit <- glm(
      Y_ij_star ~ (treatment * Y_ij) + (treatment * x1) + (treatment * x2) + (treatment *
                                                                                x3) + (x4),
      data = validation,
      family = "binomial"
    )
    
    ## Fitting P(0, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 0, treatment = 1)
    p_0_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 1, treatment = 1)
    p_1_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(0, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 0, treatment = 0)
    p_0_0 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4, cluster_num) |>
      mutate(Y_ij = 1, treatment = 0)
    p_1_0 <- predict(expit, newdata = temp, type = 'response')
  } else{
    p_0_1 <- sum(validation$Y_ij_star * (1 - validation$Y_ij) * validation$treatment) / sum((1 - validation$Y_ij) * validation$treatment)
    
    p_1_1 <- sum(validation$Y_ij_star * validation$Y_ij * validation$treatment) / sum(validation$Y_ij * validation$treatment)
    
    p_0_0 <- sum(validation$Y_ij_star * (1 - validation$Y_ij) * (1 - validation$treatment)) / sum((1 - validation$Y_ij) * (1 - validation$treatment))
    
    p_1_0 <- sum(validation$Y_ij_star * validation$Y_ij * (1 - validation$treatment)) / sum(validation$Y_ij * (1 - validation$treatment))
  }
  
  realistic_data <- cbind(data, p_0_1, p_1_1, p_0_0, p_1_0)
  
  # Calculating the actual SSW ATE
  
  E_1 <- sum(((realistic_data$treatment * realistic_data$Y_ij_star) - (pi_hat * realistic_data$p_0_1)
  ) / (pi_hat * (
    realistic_data$p_1_1 - realistic_data$p_0_1
  ))) / nrow(realistic_data)
  E_0 <- sum((((1 - realistic_data$treatment) * realistic_data$Y_ij_star
  ) - ((1 - pi_hat) * realistic_data$p_0_0)) / ((1 - pi_hat) * (
    realistic_data$p_1_0 - realistic_data$p_0_0))) / nrow(realistic_data)
  return(E_1 - E_0)
}

boot_ATE <- function(data, parametric){
  sample_ATE <- replicate(200, {
    # Should this be by individual or by clusters? Seemingly Clusters
    sample_idx <- sample(1:length(unique(data$cluster_num)), size = 30 ,replace = T)
    boot_data <- lapply(seq_along(sample_idx), function(idx){
      data[data$cluster_num == sample_idx[idx], ] |>
        mutate(cluster_num = idx)
    })
    boot_data <- data.table::rbindlist(boot_data)
    calculate_ATE(boot_data, parametric)
  })
  return(c(quantile(sample_ATE, probs = 0.025), quantile(sample_ATE, probs = 0.975)))
}

ATE_sim_one <- function(cluster_range, parametric, ICC, independent){
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
      treatment <- rbern(1)
      cluster_num <- j
      cluster_level_data <- rbind(
        cluster_level_data,
        Y_ij_0,
        Y_ij_1,
        Y_ij_star_0,
        Y_ij_star_1,
        V_ij_0,
        V_ij_1,
        treatment,
        cluster_num
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
        V_ij_0 = if_else(treatment == 0, V_ij_0, NA),
        V_ij_1 = if_else(treatment == 1, V_ij_1, NA),
        Y_ij_star_0 = if_else(treatment == 0, Y_ij_star_0, NA),
        Y_ij_star_1 = if_else(treatment == 1, Y_ij_star_1, NA),
        Y_ij_star = coalesce(Y_ij_star_0, Y_ij_star_1),
        Y_ij_0 = if_else(V_ij_0 == 1, Y_ij_0, NA),
        Y_ij_1 = if_else(V_ij_1 == 1, Y_ij_1, NA),
        Y_ij = coalesce(Y_ij_0, Y_ij_1)
      )

    ATE_est[i] <- calculate_ATE(realistic_data, parametric)
    CI_boot <- boot_ATE(realistic_data, parametric)
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


# Running the factorial model to evaluate the ATE bias

small <- c(100, 300)
large <- c(500, 1000)
ICC_vals <- c(0.01, 0.1)
size_range <- list(small, large)

parameters <- expand.grid(
  size_idx = 1:2,
  para = 0:1,
  ICC = ICC_vals,
  ind = 0:1,
  stringsAsFactors = FALSE
)

## Setting up Parallel coding
RNGkind("L'Ecuyer-CMRG")
set.seed(999)
n_jobs <- nrow(parameters)

system.time(results <- mclapply(1:n_jobs, function(i){
  param <- parameters[i, ]
  ATE_sim_one(size_range[[param$size_idx]], param$para, param$ICC, param$ind)
}, mc.cores = 50)
)

# Evaluating results

results <- do.call(rbind, lapply(results, as.data.frame))
write.csv(results, file = "MC-Random-Results.csv")
