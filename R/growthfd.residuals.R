#' Compute residuals
#'
#' This function computes residuals between measured stature data 
#' and data generated from the growth model.
#'
#' @param x Vector with input ages
#' @param y Vector with target height measurements
#' @param par Parameters of the model
#' @param model FPCA growth model  
#' @return A vector of residuals
growthfd.residuals <- function(x, y, par, model) {
  p <- growthfd.evaluate(x, par, model);
  residuals <- c(p - y, par);
  return(residuals);
}
