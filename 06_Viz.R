# File Purpose: Our current reporting of values hides a lot of important details about the distribution of what we would achieve, This file is designed 
# for data viz that would be crucial for presentation, allowing to show these visualizations

library(tidyverse)
library(MASS)
library(boot)
library(extraDistr)
library(parallel)
set.seed(123)

#################### Universal Parameterers/Functions ##############################
num_clusters <- 30
num_iter <- 500

calculate_ATE_prop <- function(data){
  validation <- data |>
    filter(!is.na(Y_ij_0) | !is.na(Y_ij_1))
  # Fitting a propensity score
  ## Model assumes perfectly IID observations
  pi_hat_model <- glm(T_ij ~ x1 + x2 + x3 + x4, data = data, family = "binomial")
  pi_hat <- predict(pi_hat_model, newdata = data, type = 'response')
  
  # Fitting Probability models, of the form P(Y_ij, T_ij, X_ij)
  
  ## Removing previously computed probabilities
  data <- data |>
    dplyr::select(-starts_with("p"))
  
    expit <- glm(
      Y_ij_star ~ (T_ij * Y_ij) + (T_ij * x1) + (T_ij * x2) + (T_ij *x3) +  x4,
      data = validation,
      family = "binomial"
    )
    
    ## Fitting P(0, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4) |>
      mutate(Y_ij = 0, T_ij = 1)
    p_0_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 1, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4) |>
      mutate(Y_ij = 1, T_ij = 1)
    p_1_1 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(0, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4) |>
      mutate(Y_ij = 0, T_ij = 0)
    p_0_0 <- predict(expit, newdata = temp, type = 'response')
    ## Fitting P(1, 0, X_ij)
    temp <- data |>
      dplyr::select(x1, x2, x3, x4) |>
      mutate(Y_ij = 1, T_ij = 0)
    p_1_0 <- predict(expit, newdata = temp, type = 'response')

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

calculate_ATE_SSO <- function(data){
  pi_hat_model <- glm(T_ij ~ x1 + x2 + x3 + x4, data = data, family = "binomial")
  pi_hat <- predict(pi_hat_model, newdata = data, type = 'response')
  
  E_1 <- mean((data$T_ij * data$Y_ij_star) / pi_hat)
  E_0 <- mean(((1 - data$T_ij) * data$Y_ij_star) / (1 - pi_hat))
  return(E_1 - E_0)
}
#################### Universal Parameterers/Functions ##############################

#################### IID ##############################

ATE <- rep(NA, num_iter)
SSW <- rep(NA, num_iter)
SSO <- rep(NA, num_iter)
for(i in 1:num_iter) {
  data <- list()
  for (j in 1:num_clusters) {
    # Generating our covariate matrix X_i per cluster
    cluster_size <- rdunif(1, cluster_range[1], cluster_range[2])
    x1 <- rnorm(cluster_size, mean = 1, sd = 1)
    x2 <- rnorm(cluster_size, mean = 0.5, sd = 1)
    x3 <- rbinom(cluster_size, 1, 0.55)
    x4 <- runif(cluster_size)
    cluster_level_data <- rbind(x1, x2, x3, x4) #X_i
    # Generating our observed outcome, Y_ij, Y*_ij, and V_ij
    # Starting with Y_ij(0) and Y_ij(1)
    Y_ij_0_pi <- inv.logit((-1) + (c(0.15, 0.2, 0.15, -0.15) %*% cluster_level_data)) # Bernoulli parameter for calculating Y_ij(0)
    Y_ij_0 <- rbern(length(Y_ij_0_pi), Y_ij_0_pi)
    Y_ij_1_pi <- inv.logit((-0.25) + (c(0.15, 0.2, 0.15, -0.15) %*% cluster_level_data)) # Bernoulli parameter for calculating Y_ij(1)
    Y_ij_1 <- rbern(length(Y_ij_1_pi), Y_ij_1_pi)
    # Now moving onto Y*_ij(0) and Y*_ij(1)
      Y_ij_0_star_pi <-  inv.logit(-1.75 + (c(0.25, -0.25, -0.15, -0.1) %*% cluster_level_data) + (1.5 * Y_ij_0))# Bernoulli parameter for calculating Y*_ij(0)
      Y_ij_star_0 <- rbern(length(Y_ij_0_star_pi), Y_ij_0_star_pi)
      Y_ij_1_star_pi <-  inv.logit(-1.25 + (c(-0.25, -0.15, -0.25, -0.1) %*% cluster_level_data) + (2.5 * Y_ij_1)) # Bernoulli parameter for calculating Y*_ij(1)
      Y_ij_star_1 <- rbern(length(Y_ij_1_star_pi), Y_ij_1_star_pi)
    # Now, we calculate V_ij
    V_ij_0_pi <- inv.logit((-0.25) + (c(-0.5, -0.5, 0.25, -0.25) %*% cluster_level_data) + (-0.15 * Y_ij_0)) # Bernoulli parameter for calculating V_ij(0)
    V_ij_0 <- rbern(length(V_ij_0_pi), V_ij_0_pi)
    V_ij_1_pi <- inv.logit((-0.5) + (c(-0.5, -0.5, 0.25, -0.25) %*% cluster_level_data) + (0.15 * Y_ij_1))  # Bernoulli parameter for calculating V_ij(1)
    V_ij_1 <- rbern(length(V_ij_1_pi), V_ij_1_pi)
    # Finally, We randomly assign our treatment (READ: this is the biggest difference between the first file and this one
    # assume IID for now, even if this logically doesn't make sense in a cluster based setting)
    T_ij_pi <- inv.logit(0.5 + (c(0.15, 0.25, -0.1, 0.1) %*% cluster_level_data)) 
    T_ij <- rbern(length(T_ij_pi), T_ij_pi)
    cluster_level_data <- rbind(
      cluster_level_data,
      Y_ij_0,
      Y_ij_1,
      Y_ij_star_0,
      Y_ij_star_1,
      V_ij_0,
      V_ij_1,
      T_ij
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
  
  SSW[i] <- calculate_ATE_prop(realistic_data)
  SSO[i] <- calculate_ATE_SSO(realistic_data)
}

SSW <- cbind(ATE = SSW, type = "SSW")
SSO <- cbind(ATE = SSO, type = "SSO")

df <- as.data.frame(rbind(SSW, SSO))
df$ATE <- as.numeric(df$ATE)

plot_IID <- df |>
  ggplot(aes(x = type, y = ATE, color = type)) +
  geom_boxplot() +
  geom_hline(yintercept = mean(ATE), linetype = "dashed", alpha = 0.5) +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom") +
  labs(x = NULL,
       y = "Estimated ATE",
       color = NULL) 

ggsave("./images/plot_IID_poster.png", plot = plot_IID, width = 12, height = 5, dpi = 300, units = "in")
ggsave("./images/plot_IID_slide.png", plot = plot_IID, width = 6, height = 4, dpi = 300, units = "in")


#################### IID ##############################

#################### Clustered, Individual ##############################

ATE <- rep(NA, num_iter)
SSW <- rep(NA, num_iter)
SSO <- rep(NA, num_iter)
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
    b_tij <- rnorm(1, 0, sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
    T_ij_pi <- inv.logit(0.5 + (c(0.1, -0.15, -0.1, 0.2) %*% cluster_level_data) + b_tij)
    T_ij <- rbern(length(T_ij_pi), T_ij_pi)
    cluster_num <- j
    cluster_data <- rbind(
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
    data[[j]] <- t(cluster_data)
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
  
  SSW[i] <- calculate_ATE_prop(realistic_data)
  SSO[i] <- calculate_ATE_SSO(realistic_data)
}

SSW <- cbind(ATE = SSW, type = "SSW")
SSO <- cbind(ATE = SSO, type = "SSO")

df <- as.data.frame(rbind(SSW, SSO))
df$ATE <- as.numeric(df$ATE)

plot_clustered <- df |>
  ggplot(aes(x = type, y = ATE, color = type)) +
  geom_boxplot() +
  geom_hline(yintercept = mean(ATE), linetype = "dashed", alpha = 0.5) +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom") +
  labs(x = NULL,
       y = "Estimated ATE",
       color = NULL) 

ggsave("./images/plot_clustered_poster.png", plot = plot_clustered, width = 12, height = 5, dpi = 300, units = "in")
ggsave("./images/plot_clustered_slide.png", plot = plot_clustered, width = 6, height = 4, dpi = 300, units = "in")

#################### Clustered, Individual ##############################

