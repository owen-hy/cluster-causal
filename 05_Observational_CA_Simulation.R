library(tidyverse)
library(MASS)
library(boot)
library(extraDistr)
library(parallel)
library(grf)

################### FUNCTIONS #######################

calculate_ATE_prop <- function(data, parametric, pi_hat){
  validation <- data |>
    filter(!is.na(Y_ij_0) | !is.na(Y_ij_1))
  
  # Fitting Probability models, of the form P(Y_ij, T_ij, X_ij)
  
  ## Removing previously computed probabilities
  data <- data |>
    dplyr::select(-starts_with("p"))
  
  if (parametric) {
    expit <- glm(
      Y_ij_star ~ (T_ij * Y_ij) + (T_ij * x1) + (T_ij * x2) + (T_ij * x3),
      data = validation,
      family = "binomial"
    )
    
    ## Fitting P(0, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3) |>
      mutate(Y_ij = 0, T_ij = 1)
    p_0_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3) |>
      mutate(Y_ij = 1, T_ij = 1)
    p_1_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(0, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3) |>
      mutate(Y_ij = 0, T_ij = 0)
    p_0_0 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3) |>
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

boot_ATE <- function(data, parametric){
  sample_ATE <- replicate(200, {
    # Should this be by individual or by clusters? Seemingly Clusters
    sample_idx <- sample(1:length(unique(data$cluster_assign)), size = 30, replace = T)
    boot_data <- lapply(seq_along(sample_idx), function(idx){
      data[data$cluster_assign == sample_idx[idx], ] |>
        mutate(cluster_assign = idx)
    })
    boot_data <- data.table::rbindlist(boot_data)
    agg <- boot_data |>
      group_by(cluster_assign) |>
      summarize(h1 = mean(x1),
              h2 = mean(x2),
              h3 = mean(x3),
              w1 = mean(w1),
              w2 = mean(w2),
              treatment = mean(T_ij))
    propensity <- probability_forest(X = agg[, c("w1", "w2", "h1", "h2", "h3")], Y = as.factor(agg$treatment))
    pi_hat <- predict(propensity, agg[, !colnames(agg) %in% c("treatment", "cluster_assign")], type = "response")$predictions[, 2]
    calculate_ATE_prop(boot_data, parametric, pi_hat[boot_data$cluster_assign])
  })
  return(c(quantile(sample_ATE, probs = 0.025), quantile(sample_ATE, probs = 0.975)))
}

############################################

n <- 100000 # Number of People
num_cluster <- 100 # Number of Clusters
num_iter <- 150
ICC <- 0.01
parametric <- T
independent <- F


ATE <- rep(NA, num_iter)
ATE_cluster <- rep(NA, num_iter)
ATE_agg <- rep(NA, num_iter)
coverage <- rep(NA, num_iter)

for(i in 1:num_iter){
  # Generating individual + cluster level covariates
  x1 <- rnorm(n)
  x2 <- runif(n)
  x3 <- rbinom(n, 1, 0.45)
  w1 <- runif(num_cluster)
  w2 <- rbinom(num_cluster, 1, 0.55)
  
  X_ij <- cbind(x1, x2, x3) # Individual Level
  W_i <- cbind(w1, w2) # Cluster Level
  
  # Assigning clusters based upon probabilities
  
  cluster_prob <- list()
  for(j in 1:num_cluster){
    cluster_prob[[j]] <- exp((0.2 * (w1[j] + w2[j]) * (1 + x1 + x2 + x3)))
  }
  cluster_prob <- prop.table(do.call(cbind, cluster_prob), margin = 1)
  cluster_assign <- apply(cluster_prob, 1, function(p) sample(seq_along(p), size = 1, prob = p))
  
  # Aggregate summary per cluster assignment
  agg <- cbind(as.data.frame(X_ij), cluster_assign) |>
    group_by(cluster_assign) |>
    summarize(h1 = mean(x1),
              h2 = mean(x2),
              h3 = mean(x3))
  
  W_i <- cbind(W_i, agg) |>
    dplyr::select(-cluster_assign)
  
  # Assigning T_ij, Y(*)_ij, V_ij
  b <- rnorm(num_cluster, mean = 0, sd = sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
  b_ind <- b[cluster_assign]
  Y_ij_0_pi <- inv.logit((-1) + (c(0.15, 0.2, 0.15) %*% t(X_ij)) + b_ind) # Bernoulli parameter for calculating Y_ij(0)
  Y_ij_0 <- rbern(length(Y_ij_0_pi), Y_ij_0_pi)
  Y_ij_1_pi <- inv.logit((-0.25) + (c(0.15, 0.2, 0.15) %*% t(X_ij)) + b_ind) # Bernoulli parameter for calculating Y_ij(1)
  Y_ij_1 <- rbern(length(Y_ij_1_pi), Y_ij_1_pi)
    
  if(independent){
    Y_ij_0_star_pi <-  inv.logit(-1.25 + (1.5 * Y_ij_0))# Bernoulli parameter for calculating Y*_ij(0)
    Y_ij_star_0 <- rbern(length(Y_ij_0_star_pi), Y_ij_0_star_pi)
    Y_ij_1_star_pi <-  inv.logit(-1 + (2.5 * Y_ij_1)) # Bernoulli parameter for calculating Y*_ij(1)
    Y_ij_star_1 <- rbern(length(Y_ij_1_star_pi), Y_ij_1_star_pi)
  } else{
    Y_ij_0_star_pi <-  inv.logit(-1.75 + (c(0.25, -0.25, -0.15) %*% t(X_ij)) + (1.5 * Y_ij_0))# Bernoulli parameter for calculating Y*_ij(0)
    Y_ij_star_0 <- rbern(length(Y_ij_0_star_pi), Y_ij_0_star_pi)
    Y_ij_1_star_pi <-  inv.logit(-1.25 + (c(-0.25, -0.15, -0.25) %*% t(X_ij)) + (2.5 * Y_ij_1)) # Bernoulli parameter for calculating Y*_ij(1)
    Y_ij_star_1 <- rbern(length(Y_ij_1_star_pi), Y_ij_1_star_pi)
  }
  
  b <- rnorm(num_cluster, mean = 0, sd = sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
  b_ind <- b[cluster_assign]
  V_ij_0_pi <- inv.logit((-0.25) + (c(-0.5, -0.5, 0.25) %*% t(X_ij)) + (-0.15 * Y_ij_0)  + b_ind) # Bernoulli parameter for calculating V_ij(0)
  V_ij_0 <- rbern(length(V_ij_0_pi), V_ij_0_pi)
  V_ij_1_pi <- inv.logit((-0.5) + (c(-0.5, -0.5, 0.25) %*% t(X_ij)) + (0.15 * Y_ij_1)  + b_ind)  # Bernoulli parameter for calculating V_ij(1)
  V_ij_1 <- rbern(length(V_ij_1_pi), V_ij_1_pi)
  
  T_ij_pi <- inv.logit(c(0.15, 0.2, 0.1, -0.175, 0.125) %*% t(W_i))
  T_ij <- rbern(length(T_ij_pi), T_ij_pi)
  
  # Combining everything 
  true_data <- as.data.frame(cbind(
    X_ij,
    w1 = w1[cluster_assign],
    w2 = w2[cluster_assign],
    Y_ij_0, Y_ij_1,
    Y_ij_star_0, Y_ij_star_1,
    V_ij_0, V_ij_1,
    T_ij = T_ij[cluster_assign],
    cluster_assign
  ))
  
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
  
  W_i <- cbind(W_i, T_ij)
  
  # Calculating Propensity Scores

  propensity_cluster <- probability_forest(X = W_i[, c("w1", "w2")], Y = as.factor(W_i$T_ij))
  propensity_agg <- probability_forest(X = W_i[, c("w1", "w2", "h1", "h2", "h3")], Y = as.factor(W_i$T_ij))
  pi_hat_cluster <- predict(propensity_cluster, W_i[, !colnames(W_i) %in% c("T_ij", "h1", "h2", "h3")], type = "response")$predictions[, 2]
  pi_hat_agg <- predict(propensity_agg, W_i[, !colnames(W_i) %in% c("T_ij")], type = "response")$predictions[, 2]
  
  # Calculations!
  ATE[i] <- mean(true_data$Y_ij_1 - true_data$Y_ij_0)
  ATE_cluster[i] <- calculate_ATE_prop(realistic_data, parametric, pi_hat_cluster[cluster_assign])
  ATE_agg[i] <- calculate_ATE_prop(realistic_data, parametric, pi_hat_agg[cluster_assign])
  
  # Coverage will only be for the aggregate, which we would expect to do better
  boot <- boot_ATE(realistic_data, parametric)
  coverage[i] <- boot[1] < ATE[i] & ATE[i] < boot[2]
}
