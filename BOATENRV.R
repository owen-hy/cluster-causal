## File Purpose: I think it could be nice to have a single function that dictates what our estimator is doing, without the mess of the Monte Carlo
## simulation code. This file is designed to showcase this estimator in a usable manner, maybe or maybe not something I can release to the public
## to have more of an impact! Good reproducible science requires well documented code!

## BOATENRV (Been pronouncing like Boat N' RV, like a fun camping trip): Binary Observational Average Treatment Effect (with) Non Random Validation (subset)
## Name is not final, but will workshop it

# While roxygen2 notations are used here, this file is not explicitly designed for R package release (yet), just for reference

#' @param df Data Frame, Matrix, or Tibble object with silver standard outcome, gold standard outcome, treatment variable, and covariates
#' @param SS String object of the name of the measurement error prone silver standard outcome
#' @param GS String object of the name of the validated gold standard outcome
#' @param treatment String object of the name of the treatment variable
#' @param non_parametric Boolean value for method of calculating silver standard probability (Default: FALSE)
#' @param propensity_model Function to use for propensity model
#' @param classification_model Function to use for classification model
#' @param propensity_args List object of arguments (with correct corresponding parameter name) to send into propensity model. Formulas must be formatted as "Treatment ~ Covariates"
#' @param classification_args List object of arguments (with correct corresponding parameter name) to send into classification model. Formulas must be formatted as "Silver Standard ~ Covariates"
#' IMPORTANT: In *_args parameter, you must refer to the design matrix as df, as that is how it will be referred to within the function
#'
#' @examples
#' 

BOATENRV <- function(df, SS, GS, treatment,
                     non_parametric = FALSE, 
                     propensity_model, classification_model, 
                     propensity_args, classification_args) {
  # TODO: ADD ... = NULL 
  # TODO: ADD rm(list = c(...))
  
  # Creating Validation subset
  validation <- df |>
    dplyr::filter(!is.na(!!rlang::sym(GS)))
  
  # Fitting classification model
  
  if(non_parametric){
    p_0_1 <- sum(df[, SS] * (1 - df[, GS]) * df[, treatment]) / sum((1 - df[, GS]) * df[, treatment])
    p_1_1 <- sum(df[, SS] * df[, GS] * df[, treatment]) / sum(df[, GS] * df[, treatment])
    p_0_0 <- sum(df[, SS] * (1 - df[, GS]) * (1 - df[, treatment])) / sum((1 - df[, GS]) * (1 - df[, treatment]))
    p_1_0 <- sum(df[, SS] * df[, GS] * (1 - df[, treatment])) / sum(df[, GS] * (1 - df[, treatment]))
  } else{
    
  }
    
  
  
}
  
  
  
  
  
