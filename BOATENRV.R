## File Purpose: I think it could be nice to have a single function that dictates what our estimator is doing, without the mess of the Monte Carlo
## simulation code. This file is designed to showcase this estimator in a usable manner, maybe or maybe not something I can release to the public
## to have more of an impact! Good reproducible science requires well documented code!

## BOATENRV (Been pronouncing like Boat N' RV, like a fun camping trip): Binary Observational Average Treatment Effect (with) Non Random Validation (subset)
## Name is not final, but will workshop it

# While roxygen2 notations are used here, this file is not explicitly designed for R package release (yet), just for reference

#' @param SS Vector of error prone silver standard outcomes 
#' @param GS Vector of validated gold standard outcomes (Observations without validation should be NA)
#' @param treatment Vector of binary treatment assignment
#' @param propensity_X Covariates of interest in propensity model (Matrix, Tibble, or Data Frame)
#' @param classification_X Covariates of interest in classification model (Matrix, Tibble, or Data Frame)
#' @param propensity_model Function to use for propensity model
#' @param classification_model Function to use for classification model
#' @param non_parametric Boolean value for method of calculating silver standard probability (Default: FALSE)
#'
#' @examples
#' 
#' 
#' 

BOATENRV <- function(SS, GS, treatment, propensity_X, classification_X,
                     propensity_model = c("glm", "glmer", "rf"), classification_model = c("glm", "glmer", "rf"), 
                     non_parametric = FALSE) {
  # TODO: ADD ... = NULL 
  # TODO: ADD rm(list = c(...))
  
  # Arguments for Model, matching it correspondingly
  propensity_model <- match.arg(propensity_model, c("glm", "glmer", "rf"))
  classification_model <- match.arg(classification_model, c("glm", "glmer", "rf"))
  prop_formula <- T
  if(popensity_model == "glm"){
    prop_model <- glm
  } else if(propensity_model == "glmer"){
    prop_model <- lme4::glmer
  } else{
    prop_model <- grf::probability_forest
    prop_formula <- F
  }
  
  class_formula <- T
  if(classification_model == "glm"){
    class_model <- glm
  } else if(classification_model == "glmer"){
    class_model <- lme4::glmer
  } else{
    class_model <- grf::probability_forest
    class_formula <- F
  }
  
  # Data Pre-Processing
  
  df <- cbind(SS, GS, treatment, propensity_X)
  
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
  validation <- df |>
    dplyr::filter(!is.na(GS))
  
  # Fitting classification model
  
  if(non_parametric){
    p_0_1 <- sum(validation$SS * (1 - validation$GS) * validation$treatment) / sum((1 - validation$GS) * validation$treatment)
    p_1_1 <- sum(validation$SS * validation$GS * validation$treatment) / sum(validation$GS * validation$treatment)
    p_0_0 <- sum(validation$SS * (1 - validation$GS) * (1 - validation$treatment)) / sum((1 - validation$GS) * (1 - validation$treatment))
    p_1_0 <- sum(validation$SS * validation$GS * (1 - validation$treatment)) / sum(validation$GS * (1 - validation$treatment))
  } else{
    class_names <- colnames(classification_X)
    if(class_formula){ # This is the glm/glmer scenario, formula can be used
      
      # TODO: PLEASE, make sure that we are including treatment, goldstandard, and interactions within this model
      
    } else{ # This is the random forest scenario, X and Y must be specified in specific way
      
      # TODO: PLEASE, make sure that we are including treatment, goldstandard, and interactions within this model
      
    }
    
    
  }
  
  # Fitting Propensity Model
  
  if(prop_formula){ # glm/glmer scenario for propensity score, formula can be used
    # TODO: fit pi_hat
  } else{ # Random forest scenario for propensity score, matrix notation
    # TODO: fit pi_hat
  }
  
  # Constructing the estimator itself
  E_1 <- sum(((treatment * SS) - (pi_hat * p_0_1)) / (pi_hat * (p_1_1 - p_0_1))) / length(SS)
    
  E_0 <- sum((((1 - treatment) * SS) - ((1 - pi_hat) * p_0_0)) / ((1 - pi_hat) * (p_1_0 - p_0_0))) / length(SS)
    
  return(E_1 - E_0)
  
}
  
  
  
  
  
