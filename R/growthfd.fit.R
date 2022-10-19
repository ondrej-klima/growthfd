#' Fit a FPCA Growth Curve Model to measurements of a single individual
#'
#' This function fits a model to the given measured data of a single individual.
#'
#' @param model FPCA growth model to be fitted
#' @param age Age at measured data points
#' @param height Height at at measured data points
#' @param nprint Verbosity
#' @return An optimization result object
#' @example man/examples/fit.R
#' @export
growthfd.fit <- function(model, age, height, nprint=1) {
  #defaultW <- getOption("warn") 
  #options(warn = -1) 
  
  npar <- sum(model$scores.elements);
  suppressWarnings(r<-minpack.lm::nls.lm(
    par=rep(0,npar), 
    fn = growthfd.residuals, 
    x = age, 
    y = height, 
    model=model, 
    control = minpack.lm::nls.lm.control(nprint=nprint),
    upper=rep(20,npar), 
    lower=rep(-20,npar)))
  
#    upper=c(10,10,10,3,3,3,10,10,10,1000,3,3), #rep(10,npar), 
#    lower=-c(10,10,10,3,3,3,10,10,10,1000,3,3))) #rep(-10,npar)))
  #options(warn = defaultW)
  return(r)
}
