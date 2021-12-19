#' Generate a Curve Function
#'
#' This function generates a growth curve function based on given model 
#' and parameters, describing the growth phase and amplitude.
#'
#' @param par Phase and amplitude parameters
#' @param model FPCA growth model 
#' @return FDA function object
#' @export
growthfd.std <- function(par, model) {
  library('fda')
  nwarpscores <- model$scores.elements[1];
  ngrowthscores <- model$scores.elements[2];
  
  warpscores <- par[1:nwarpscores];
  growthscores <- par[(nwarpscores+1):length(par)];
  
  growthfd <- model$growthfpca$meanfd;
  for(i in 1:ngrowthscores) {
    growthfd <- growthfd + growthscores[i]*sqrt(model$growthfpca$values[i])*model$growthfpca$harmonics[i];
  }
  
  warpfd <- model$warpfpca$meanfd;
  for(i in 1:nwarpscores) {
    warpfd <- warpfd + warpscores[i]*sqrt(model$warpfpca$values[i])*model$warpfpca$harmonics[i];
  }
  
  return(register.newfd(growthfd, warpfd));
}

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
  p <- growthfd.evaluate(x, par, model);
  residuals <- c(p - y, par);
  return(residuals);
}

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
  nscores <- sum(model$scores.elements);
  scores <- par[1:sum(nscores)];
  
  growthfd <- growthfd.std(scores, model);
  return(fda::eval.fd(x, growthfd, deriv));
}

#' Plot a Growth Curve
#'
#' This function plots a stature, velocity or acceleration curve.
#'
#' @param par Parameters of the model
#' @param model FPCA growth model
#' @param deriv Path to the input file 
#' @param from The lower age limit
#' @param to The upper age limit
#' @export
growthfd.plot <- function(model, par, deriv=0, from=0, to=18) {
  x <- seq(from,to,0.05);
  ylab <- switch(deriv+1, "Stature (cm)", "Velocity (cm per year)", "Acceleration (cm per year sqr.)")
  plot(x, growthfd.evaluate(x, par, model, deriv), xlab = "Age (years)", ylab = ylab, type="l", col="blue")
}

#' Fit a FPCA Growth Curve Model to the Data
#'
#' This function fits a model to the given measured data.
#'
#' @param model FPCA growth model to be fitted
#' @param age Age at measured data points
#' @param height Height at at measured data points
#' @param nprint Verbosity
#' @return An optimization result object
#' @example man/examples/fit.R
#' @export
growthfd.fit <- function(model, age, height, nprint=1) {
  npar <- sum(model$scores.elements);
  return(minpack.lm::nls.lm(par=rep(0,npar), fn = growthfd.residuals, x = age, y = height, model=model, control = minpack.lm::nls.lm.control(nprint=1), upper=rep(3,npar), lower=rep(-3,npar)))
}

growthfd.growthfd <- function(data, x, y, id, model, verbose=1) {
  
}