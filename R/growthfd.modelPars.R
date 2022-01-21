#' Standardized model scores
#' 
#' This function returns model parameters for the individuals used for training
#' the model.
#' 
#' @param model FPCA-based growth model
#' @return Matrix containing the scores
#' @export
growthfd.modelPars <- function(model) {
  warpn <- model$scores.elements[1];
  growthn <- model$scores.elements[2];
  
  n <- nrow(model$growthfpca$scores)
  
  values <- matrix(rep(model$warpfpca$values[1:warpn], each=n), nrow=n);
  warpscores <- model$warpfpca$scores / sqrt(values); 
  
  values <- matrix(rep(model$growthfpca$values[1:growthn], each=n), nrow=n);
  growthscores <- model$growthfpca$scores / sqrt(values); 
  
  return(cbind(warpscores, growthscores));
}