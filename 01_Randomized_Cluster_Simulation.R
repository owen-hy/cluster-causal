# File Purpose: This R script is designed to replicate the result of Dane's estimator in a simulation study
# This is done under the ideal situation where all clusters are randomized (Thus no propensity model is needed)

library(tidyverse)
library(MASS)
library(boot)
library(extraDistr)
library(geepack)
set.seed(123)

start <- Sys.time()

num_iter <- 100 # Dane used 5000, just doing 100 until I learn how to parallel code
num_clusters <- 30
ATE <- rep(NA, num_iter)
ATE_est <- rep(NA, num_iter)

ATE_sim_one <- function(cluster_range, parametric, ICC, independent){
  for(i in 1:num_iter) {
    data <- list()
    for (j in 1:num_clusters) {
      # Generating our covariate matrix X_i per cluster
      cluster_size <- rdunif(1, cluster_range[1], cluster_range[2])
      x1 <- mvrnorm(mu = rep(1, cluster_size),
                    Sigma = diag(x = cluster_size))
      x2 <- mvrnorm(mu = rep(0.5, cluster_size), Sigma = ((0.5 * diag(x = cluster_size)) + (
        0.05 * matrix(
          data = 1,
          nrow = cluster_size,
          ncol = cluster_size
        )
      )))
      x3 <- rbern(cluster_size, prob = 0.55)
      x4 <- rep(runif(1), cluster_size)
      cluster_level_data <- rbind(x1, x2, x3, x4) #X_i
      # Generating our observed outcome, Y_ij, Y*_ij, and V_ij
      # Starting with Y_ij(0) and Y_ij(1)
      b_i_y <- rnorm(1, 0, sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
      Y_ij_0_pi <- inv.logit((-1) + (c(0.15, 0.2, 0.15, -0.15) %*% cluster_level_data) + b_i_y) # Bernoulli parameter for calculating Y_ij(0)
      Y_ij_0 <- sapply(Y_ij_0_pi, function(p)
        rbern(1, p))
      Y_ij_1_pi <- inv.logit((-0.25) + (c(0.15, 0.2, 0.15, -0.15) %*% cluster_level_data) + b_i_y) # Bernoulli parameter for calculating Y_ij(1)
      Y_ij_1 <- sapply(Y_ij_1_pi, function(p)
        rbern(1, p))
      # Now moving onto Y*_ij(0) and Y*_ij(1)
      b_i_ystar <- 0 # READ: Paper stated that they kept this at zero?
      if(independent){
        Y_ij_0_star_pi <-  inv.logit(-1.25 + (1.5 * Y_ij_0) + b_i_ystar)# Bernoulli parameter for calculating Y*_ij(0)
        Y_ij_star_0 <- sapply(Y_ij_0_star_pi, function(p)
          rbern(1, p))
        Y_ij_1_star_pi <-  inv.logit(-1 + (2.5 * Y_ij_1) + b_i_ystar) # Bernoulli parameter for calculating Y*_ij(1)
        Y_ij_star_1 <- sapply(Y_ij_1_star_pi, function(p)
          rbern(1, p))
      } else{
        Y_ij_0_star_pi <-  inv.logit(-1.75 + (c(0.25, -0.25, -0.15, -0.1) %*% cluster_level_data) + (1.5 * Y_ij_0) + b_i_ystar)# Bernoulli parameter for calculating Y*_ij(0)
        Y_ij_star_0 <- sapply(Y_ij_0_star_pi, function(p)
          rbern(1, p))
        Y_ij_1_star_pi <-  inv.logit(-1.25 + (c(-0.25, -0.15, -0.25, -0.1) %*% cluster_level_data) + (2.5 * Y_ij_0) + b_i_ystar) # Bernoulli parameter for calculating Y*_ij(1)
        Y_ij_star_1 <- sapply(Y_ij_1_star_pi, function(p)
          rbern(1, p))
      }
      # Now, we calculate V_ij
      b_iv <- rnorm(1, 0, sqrt((ICC * (pi^2)) / (3 * (1 - ICC))))
      V_ij_0_pi <- inv.logit((-0.25) + (c(-0.5, -0.5, 0.25, -0.25) %*% cluster_level_data) + (-0.15 * Y_ij_0)  + b_iv) # Bernoulli parameter for calculating V_ij(0)
      V_ij_0 <- sapply(V_ij_0_pi, function(p)
        rbern(1, p))
      V_ij_1_pi <- inv.logit((-0.5) + (c(-0.5, -0.5, 0.25, -0.25) %*% cluster_level_data) + (0.15 * Y_ij_1)  + b_iv)  # Bernoulli parameter for calculating V_ij(1)
      V_ij_1 <- sapply(V_ij_1_pi, function(p)
        rbern(1, p))
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

    validation <- realistic_data |>
      filter(!is.na(Y_ij_0) | !is.na(Y_ij_1))
    pi_hat <- mean(realistic_data$treatment)

    # Fitting logistic model
    expit <- geeglm(
      Y_ij_star ~ (treatment * Y_ij) + (treatment * x1) + (treatment * x2) + (treatment *
                                                                                x3) + (x4),
      data = validation,
      family = "binomial",
      id = cluster_num
    )

    # Fitting Probability models, of the form P(Y_ij, T_ij, X_ij)

    if (parametric) {
      ## Fitting P(0, 1, X_ij)
      temp <- realistic_data |>
        dplyr::select(x1, x2, x3, x4, cluster_num) |>
        mutate(Y_ij = 0, treatment = 1)
      p_0_1 <- predict(expit, newdata = temp, type = 'response')
      ## Fitting P(1, 1, X_ij)
      temp <- realistic_data |>
        dplyr::select(x1, x2, x3, x4, cluster_num) |>
        mutate(Y_ij = 1, treatment = 1)
      p_1_1 <- predict(expit, newdata = temp, type = 'response')
      ## Fitting P(0, 0, X_ij)
      temp <- realistic_data |>
        dplyr::select(x1, x2, x3, x4, cluster_num) |>
        mutate(Y_ij = 0, treatment = 0)
      p_0_0 <- predict(expit, newdata = temp, type = 'response')
      ## Fitting P(1, 0, X_ij)
      temp <- realistic_data |>
        dplyr::select(x1, x2, x3, x4, cluster_num) |>
        mutate(Y_ij = 1, treatment = 0)
      p_1_0 <- predict(expit, newdata = temp, type = 'response')
    } else{
      p_0_1 <- sum(validation$Y_ij_star * (1 - validation$Y_ij) * validation$treatment) / sum((1 - validation$Y_ij) * validation$treatment)

      p_1_1 <- sum(validation$Y_ij_star * validation$Y_ij * validation$treatment) / sum(validation$Y_ij * validation$treatment)

      p_0_0 <- sum(validation$Y_ij_star * (1 - validation$Y_ij) * (1 - validation$treatment)) / sum((1 - validation$Y_ij) * (1 - validation$treatment))

      p_1_0 <- sum(validation$Y_ij_star * validation$Y_ij * (1 - validation$treatment)) / sum(validation$Y_ij * (1 - validation$treatment))
    }

    realistic_data <- cbind(realistic_data, p_0_1, p_1_1, p_0_0, p_1_0)

    # Calculating the actual SSW ATE

    E_1 <- sum(((realistic_data$treatment * realistic_data$Y_ij_star) - (pi_hat * realistic_data$p_0_1)
    ) / (pi_hat * (
      realistic_data$p_1_1 - realistic_data$p_0_1
    ))) / nrow(realistic_data)
    E_0 <- sum((((1 - realistic_data$treatment) * realistic_data$Y_ij_star
    ) - ((1 - pi_hat) * realistic_data$p_0_0)) / ((1 - pi_hat) * (
      realistic_data$p_1_0 - realistic_data$p_0_0
    ))) / nrow(realistic_data)
    ATE_est[i] <- E_1 - E_0
  }
}

end <- Sys.time()

end - start


# Running the factorial model to evaluate the ATE bias

small <- c(100, 300)
large <- c(500, 1000)
size <- list(small, large)



# For troubleshooting
ggplot(realistic_data) +
  geom_density(aes(x = p_0_0), color = "red") +
  geom_density(aes(x = p_0_1), color = "blue") +
  geom_density(aes(x = p_1_1), color = "green") +
  geom_density(aes(x = p_1_0), color = "black") +
  labs(title = "Density of Probabilities",
       subtitle = "Denominators are P(1, 1) (Green) - P(0, 1) (Blue) & P(1, 0) (Black)- P(0, 0) (Red)",
       x = "Probability") +
  theme_bw()
