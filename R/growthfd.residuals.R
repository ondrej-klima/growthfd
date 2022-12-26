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
#' @export
growthfd.residuals <- function(x, y, par, model) {
  suppressWarnings(p <- growthfd.evaluate(x, par, model));
  #d <- diff(growthfd.evaluate(seq(7.5, 20, 0.5), par, model))
  #d[d > 0] <- 0
  # residuals <- c(10 * (p - y), par, 100 * d * sum(d < 0));
  #residuals <- c(10 * (p - y), par, 10* d);
  residuals <- c(10 * (p - y), par);
  return(residuals);
}
