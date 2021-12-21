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
  
  return(fda::register.newfd(growthfd, warpfd));
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
  npar <- sum(model$scores.elements);
  return(minpack.lm::nls.lm(par=rep(0,npar), fn = growthfd.residuals, x = age, y = height, model=model, control = minpack.lm::nls.lm.control(nprint=1), upper=rep(3,npar), lower=rep(-3,npar)))
}

#' Fit a FPCA Growth Curve Model to a population
#'
#' This function fits a model to the given measured data of a population.
#'
#' @param model FPCA growth model to be fitted
#' @param data Data frame containing age, height and id of individuals
#' @param x Age at measured data points
#' @param y Height at measured data points
#' @param id Corresponding individual's id at measured data points
#' @param verbose Verbosity
#' @return List containing individuals id and model 
#' @example man/examples/main.R
#' @export
growthfd <- function(data, x, y, id, model, verbose=1) {
  mcall <- match.call()
  x.na <- as.numeric(eval(mcall$x, data))
  y.na <- as.numeric(eval(mcall$y, data))
  id.na <- as.factor(eval(mcall$id, data))
    
  msk <- !is.na(x.na) & !is.na(y.na)
  x <- x.na[msk]
  y <- y.na[msk]
  id <- id.na[msk]
  
  msk <- !is.na(x.na)
  x.na <- x.na[msk]
  y.na <- y.na[msk]
  id.na <- id.na[msk]

  ids <- levels(id)
  n <- length(ids)
  
  scores <- matrix(NA, n, sum(model$scores.elements))
  milestones <- matrix(NA, n, 6)
  colnames(milestones) <- c('apv', 'vpv', 'hpv', 'atf', 'vtf', 'htf')
  m_x <- seq(7, 18, 0.05)
  fitted <- data.frame('id'=factor(), 'fitted'=double(), 'residuals'=double())
  
  sampling <- seq(0, 18, 0.25)
  m <- length(sampling)
  stature <- matrix(NA, n, m)
  velocity <- matrix(NA, n, m)
  acceleration <- matrix(NA, n, m)
  
  
  for(i in seq_along(ids)) {
    msk <- id == ids[i]
    
    if(verbose > 0) {
      message(sprintf("Processing individual with id '%s' (%d/%d), containing %d measurements\n", ids[i], i, n, sum(msk)))
    }
    
    fit <- growthfd.fit(model, x[msk], y[msk], verbose)
    scores[i,] <- fit$par
    
    # Computation of growth milestones
    m_y <- growthfd.evaluate(m_x, fit$par, model, deriv=1)
    xyi <- data.frame('x'=m_x, 'y'=m_y)
    
    peak.xyi <- sitar::getPeak(xyi)
    milestones[i, 'apv'] <- peak.xyi[1]
    milestones[i, 'vpv'] <- peak.xyi[2]
    takeoff.xyi <- sitar::getTakeoff(xyi)
    milestones[i, 'atf'] <- takeoff.xyi[1]
    milestones[i, 'vtf'] <- takeoff.xyi[2]
    
    if(!is.na(milestones[i, 'apv'])) {
      milestones[i, 'hpv'] <- growthfd.evaluate(milestones[i, 'apv'], fit$par, model)
    }

    if(!is.na(milestones[i, 'atf'])) {
      milestones[i, 'htf'] <- growthfd.evaluate(milestones[i, 'atf'], fit$par, model)
    }
    
    # Evaluation of fits and residuals
    msk.na <- id.na == ids[i]
    
    f <- growthfd.evaluate(x.na[msk.na], fit$par, model)
    r <- f - y[msk.na]
    fitted <- rbind(fitted, data.frame('id'=id[msk.na], 'fitted'=f, 'residuals'=r))
    
    # Evaluations of stature, velocity and acceleration
    stature[i,] <- growthfd.evaluate(sampling, fit$par, model)
    velocity[i,] <- growthfd.evaluate(sampling, fit$par, model)
    acceleration[i,] <- growthfd.evaluate(sampling, fit$par, model)
  }
  
  colnames(fitted) <- c('id', 'fitted', 'residuum')
  return(list('ids' = ids, 'scores' = scores, 'milestones' = milestones, 'fitted' = fitted, 'stature' = stature, 'velocity' = velocity, 'acceleration' = acceleration))
}