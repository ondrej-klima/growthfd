#' Generate a Discrete Growth Curve
#'
#' This function evaluates a curve function for given ages. Depending 
#' on a degree of derivation, the function produces stature, velocity 
#' or acceleration curve.
#' 
#' @param x Ages to be evaluated
#' @param par Parameters of the model
#' @param model FPCA growth model
#' @param deriv Path to the input file 
#' @return Y-values of the evaluated curve
#' @export
growthfd.evaluate <- function(x, par, model, deriv=0) {
  defaultW <- getOption("warn") 
  options(warn = -1) 
  
  nscores <- sum(model$scores.elements);
  scores <- par[1:sum(nscores)];
  
  growthfd <- growthfd.std(scores, model);
  r <- fda::eval.fd(x, growthfd, deriv)
  options(warn = defaultW)
  return(r);
}
