## File Purpose: I think it could be nice to have a single function that dictates what our estimator is doing, without the mess of the Monte Carlo
## simulation code. This file is designed to showcase this estimator in a usable manner, maybe or maybe not something I can release to the public
## to have more of an impact! Good reproducible science requires well documented code!

## BOATENRV (Been pronouncing like Boat N' RV, like a fun camping trip): Binary Observational Average Treatment Effect (with) Non Random Validation (subset)
## Name is not final, but will workshop it

# While roxygen2 notations are used here, this file is not explicitly designed for R package release (yet), just for reference

#' Average Treatment Effect (ATE) estimator for binary outcome observational trials with measurement error and non-random validation subset
#' 
#' @description
#' The function calculates an ATE from data where an error prone (silver standard) and validated subset (gold standard) binary outcome exists.
#' 
#' @param SS Vector of error prone silver standard outcomes 
#' @param GS Vector of validated gold standard outcomes (Observations without validation should be NA)
#' @param treatment Vector of binary treatment assignment
#' @param propensity_X Covariates of interest in propensity model (Matrix, Tibble, or Data Frame)
#' @param classification_X Covariates of interest in classification model (Matrix, Tibble, or Data Frame)
#' @param propensity_model Function to use for propensity model
#' @param classification_model Function to use for classification model
#' @param non_parametric Boolean value for method of calculating silver standard probability (Default: FALSE)
#' @param propensity_formula Optional argument for propensity model formula as a string, for more complex predictors such as interactions. Default formula will be Treatment ~ Covariates in propensity_X. Only applicable for propensity_model == "glm" and propensity_model == "glmer"
#' @param classification_formula Optional argument for classification model formula as a string, for more complex predictors such as interactions. Default formula will be SS ~ treatment * (GS + classification_X). Only applicable for propensity_model == "glm" and propensity_model == "glmer"
#' @param cluster Vector of cluster assignment. If left empty, IID will be assumed. Cluster is required when model == "glmer" and optional when model == "rf"
#' 
#' @return List with estimated ATE and formulas for both propensity and classification model
#' 
#' @examples
#' 
#' 
#' 

BOATENRV <- function(SS, GS, treatment, propensity_X, classification_X,
                     propensity_model = c("glm", "glmer", "rf"), classification_model = c("glm", "glmer"), 
                     non_parametric = FALSE, propensity_formula = NULL, classification_formula = NULL, cluster = NULL) {
  # TODO: ADD ... = NULL 
  # TODO: ADD rm(list = c(...))
  
  # Arguments for Model, matching it correspondingly
  propensity_model <- match.arg(propensity_model, c("glm", "glmer", "rf"))
  classification_model <- match.arg(classification_model, c("glm", "glmer"))
  i_p_formula <- T
  if(propensity_model == "glm"){
    prop_model <- glm
  } else if(propensity_model == "glmer"){
    if(is.null(cluster)){
      warning("Cluster argument not provided, defaulting to glm")
      prop_model <- glm
    } else{
      prop_model <- lme4::glmer 
    }
  } else{
    prop_model <- grf::probability_forest
    i_p_formula <- F
  }
  
  if(classification_model == "glm"){
    class_model <- glm
  } else if(classification_model == "glmer"){
    if(is.null(cluster)){
      warning("Cluster argument not provided, defaulting to glm")
      class_model <- glm
      classification_model <- "glm"
    } else{
      class_model <- lme4::glmer 
    }
  } 
  
  # Data Pre-Processing
  df <- cbind(SS, GS, treatment, propensity_X)
  if(!is.null(cluster)){
    df <- cbind(df, cluster)
  }
  
  ## Propensity and Classification X might have overlapping covariates, accounting for that when binding
  
  if(!non_parametric){ # We would only need this merge if we are in a non_parametric scenario
    length_prop <- length(colnames(propensity_X))
    length_class <- length(colnames(classification_X))
    
    if(length_prop == 0 | length_class == 0){
      stop("Please add Column Names to Classification and/or Propensity Covariates")
    }
    
    missing <- colnames(classification_X)[!colnames(classification_X) %in% colnames(propensity_X)]
    df <- cbind(df, classification_X[, missing])
  }
  
  
  # Creating Validation subset
  validation <- as.data.frame(df) |>
    dplyr::filter(!is.na(GS))
  
  # Fitting classification model
  
  if(non_parametric){
    p_0_1 <- sum(validation$SS * (1 - validation$GS) * validation$treatment) / sum((1 - validation$GS) * validation$treatment)
    p_1_1 <- sum(validation$SS * validation$GS * validation$treatment) / sum(validation$GS * validation$treatment)
    p_0_0 <- sum(validation$SS * (1 - validation$GS) * (1 - validation$treatment)) / sum((1 - validation$GS) * (1 - validation$treatment))
    p_1_0 <- sum(validation$SS * validation$GS * (1 - validation$treatment)) / sum(validation$GS * (1 - validation$treatment))
  } else{
    class_names <- colnames(classification_X)
    if(classification_model == "glm") {
        if(is.null(classification_formula)){
          c_formula <- as.formula(paste(
            "SS ~ treatment * ", "(", paste(c("GS", colnames(propensity_X)), collapse = " + "), ")", sep = ""
          ))
        } else{
          c_formula <- as.formula(classification_formula)
        }
    } else{ # glmer scenario
        if(is.null(classification_formula)){
          c_formula <- as.formula(paste(
            "SS ~ (1 | cluster) + treatment * ",
            "(",
            paste(c("GS", colnames(propensity_X)), collapse = " + "),
            ")",
            sep = ""
          ))
        } else{
          c_formula <- as.formula(classification_formula) 
        }
    }
    classification <- class_model(c_formula, data = validation, family = "binomial")
    
    ## Fitting P(0, 1, X_ij)
    temp <- df[, colnames(propensity_X)] |>
      mutate(GS = 0, treatment = 1)
    p_0_1 <- predict(classification, newdata = temp, type = 'response')
    ## Fitting P(1, 1, X_ij)
    temp <- df[, colnames(propensity_X)] |>
      mutate(GS = 1, treatment = 1)
    p_1_1 <- predict(classification, newdata = temp, type = 'response')
    ## Fitting P(0, 0, X_ij)
    temp <- df[, colnames(propensity_X)] |>
      mutate(GS = 0, treatment = 0)
    p_0_0 <- predict(classification, newdata = temp, type = 'response')
    ## Fitting P(1, 0, X_ij)
    temp <- df[, colnames(propensity_X)] |>
      mutate(GS = 1, treatment = 0)
    p_1_0 <- predict(classification, newdata = temp, type = 'response')
  }

  
  # Fitting Propensity Model
  
  if(i_p_formula){ # glm/glmer scenario for propensity score, formula can be used
    if(propensity_model == "glm"){ 
      if(is.null(propensity_formula)){
        p_formula <- as.formula(paste("treatment ~ ", paste(colnames(propensity_X), collapse = " + ")))
      } else{
        p_formula <- as.formula(propensity_formula)
      }
    } else{
      if(is.null(propensity_formula)){
        p_formula <- as.formula(paste("treatment ~ ", paste(c("(1 | cluster)", colnames(propensity_X)), collapse = " + ")))
      } else{
        p_formula <- as.formula(propensity_formula)
      }
    }
    propensity <- prop_model(p_formula, data = df, family = "binomial")
    pi_hat <- predict(propensity, df, type = "response")
  } else{ # Random forest scenario for propensity score, matrix notation
    if(!is.null(propensity_formula)){
      warning("Propensity formula is not applicable to random forest model, please specify all covariates (including interactions) in propensity_X")
    }
    propensity <- prop_model(X = propensity_X, Y = as.factor(treatment), clusters = cluster)
    pi_hat <- predict(propensity, cbind(propensity_X, cluster))
    p_formula <- NULL
  }
  
  # Constructing the estimator itself
  E_1 <- sum(((treatment * SS) - (pi_hat * p_0_1)) / (pi_hat * (p_1_1 - p_0_1))) / length(SS)
    
  E_0 <- sum((((1 - treatment) * SS) - ((1 - pi_hat) * p_0_0)) / ((1 - pi_hat) * (p_1_0 - p_0_0))) / length(SS)
    
  return(list(ATE = E_1 - E_0,
         Propensity_Formula = p_formula,
         Classification_Formula = c_formula))
  
}
  
  
  
  
  
