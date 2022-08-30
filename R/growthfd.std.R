#' Generate a Curve Function
#'
#' This function generates a growth curve function based on given model 
#' and parameters, describing the growth phase and amplitude.
#'
#' @param par Phase and amplitude parameters
#' @param model FPCA growth model 
#' @return FDA function object
#' @import fda
#' @export
growthfd.std <- function(par, model) {
  nwarpscores <- model$scores.elements[1];
  ngrowthscores <- model$scores.elements[2];
  growthscores <- par[(nwarpscores+1):length(par)];
  
  growthfd <- model$growthfpca$meanfd;
  for(i in 1:ngrowthscores) {
    growthfd <- growthfd + growthscores[i]*sqrt(model$growthfpca$values[i])*model$growthfpca$harmonics[i];
  }
  
  return(fda::register.newfd(growthfd, growthfd.warpfd(par, model),"direct"));
}
