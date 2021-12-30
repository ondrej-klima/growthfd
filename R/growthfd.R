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
  
  return(fda::register.newfd(growthfd, growthfd.warpfd(par, model)));
}


#' @import fda
growthfd.warpfd <- function(par, model) {
  nwarpscores <- model$scores.elements[1];
  warpscores <- par[1:nwarpscores];
  
  warpfd <- model$warpfpca$meanfd;
  for(i in 1:nwarpscores) {
    warpfd <- warpfd + warpscores[i]*sqrt(model$warpfpca$values[i])*model$warpfpca$harmonics[i];
  }
  
  return(warpfd)
}

#' @import fda
growthfd.warpfdInv <- function(par, model) {
  ages <- seq(0, 18, 0.05)
  r <- eval.fd(ages, growthfd.warpfd(par, model))
  rng <- c(0,18)
  r[1] <- 0
  r[r < 0] <- 0
  r[r > 18] <- 17.999999999
  r[length(r)] <- 18
    
  norder <- 6
  nage <- length(r)
  nbasis <- nage + norder - 2
  wbasis <- create.bspline.basis(rangeval=rng, nbasis=nbasis, norder=norder, breaks=r)
  
  # starting values for coefficient
  cvec0 <- matrix(0,nbasis,1)
  Wfd0 <- fd(cvec0, wbasis)
  # set up functional parameter object
  Lfdobj <- 3 # penalize curvature of acceleration
  #lambda <- 10^(-0.5) # smoothing parameter
  lambda <- 5e-2 # smoothing parameter 1e-2 5e-2
  growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
  
  wbasis <- create.bspline.basis(rangeval=c(0,18), nbasis=nbasis, norder=norder, breaks=ages)
  Wfd0 <- fd(cvec0, wbasis)
  growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
  return(smooth.basis(ages, eval.fd(ages, smooth.basis(r, ages, growfdPar)$fd), growfdPar)$fd)
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
#' @param bounds Limitation of the interval for milestones estimation, 'negative' or 'inverse'
#' @return List containing individuals id and model 
#' @example man/examples/main.R
#' @export
growthfd <- function(data, x, y, id, model, verbose=1, bounds='negative') {
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
  fitted <- data.frame('id'=factor(), 'fitted'=double(), 'residuals'=double())
  
  sampling <- seq(0, 18, 0.25)
  m <- length(sampling)
  stature <- matrix(NA, n, m)
  velocity <- matrix(NA, n, m)
  acceleration <- matrix(NA, n, m)
  
  # Determine model apv and atf
  m_x <- seq(7, 18, 0.05)
  m_y <- growthfd.evaluate(m_x, rep(0, sum(model$scores.elements)), model, deriv=1)
  xyi <- data.frame('x'=m_x, 'y'=m_y)
  m_apv <- peak.xyi <- sitar::getPeak(xyi)[1]
  m_atf <- sitar::getTakeoff(xyi)[1]
  message(sprintf('Model apv=%f, atf=%f', m_apv, m_atf))
  
  
  for(i in seq_along(ids)) {
    msk <- id == ids[i]
    
    if(verbose > 0) {
      message(sprintf("Processing individual with id '%s' (%d/%d), containing %d measurements\n", ids[i], i, n, sum(msk)))
    }
    
    fit <- growthfd.fit(model, x[msk], y[msk], verbose)
    scores[i,] <- fit$par
    if(bounds == 'negative') {
      warpfd <- growthfd.warpfd(-fit$par, model)
    }
    else {
      warpfd <- growthfd.warpfdInv(fit$par, model)
    }
    w_m <- fda::eval.fd(c(m_apv, m_atf), warpfd)
    #w_apv <- 2*m_apv - w_m[1]
    #w_atf <- 2*m_atf - w_m[2]
    w_apv <- w_m[1]
    w_atf <- w_m[2]
    message(sprintf('Warped apv=%f, atf=%f', w_apv, w_atf))
    
    # Computation of growth milestones
    # apv
    m_x_apv <- seq(w_apv-2, w_apv+2, 0.05)
    # m_x_apv <- seq(7, 18, 0.05)
    m_y_apv <- growthfd.evaluate(m_x_apv, fit$par, model, deriv=1)
    xyi_apv <- data.frame('x'=m_x_apv, 'y'=m_y_apv)
    
    peak.xyi_apv <- sitar::getPeak(xyi_apv)
    milestones[i, 'apv'] <- peak.xyi_apv[1]
    milestones[i, 'vpv'] <- peak.xyi_apv[2]
    
    # atf
    m_x_atf <- seq(w_atf-2, w_apv+2, 0.05)
    # m_x_atf <- seq(7, 18, 0.05)
    m_y_atf <- growthfd.evaluate(m_x_atf, fit$par, model, deriv=1)
    xyi <- data.frame('x'=m_x_atf, 'y'=m_y_atf)
    
    takeoff.xyi <- sitar::getTakeoff(xyi)
    milestones[i, 'atf'] <- takeoff.xyi[1]
    milestones[i, 'vtf'] <- takeoff.xyi[2]
    
    message(sprintf('Refined apv=%f, atf=%f', milestones[i, 'apv'], milestones[i, 'atf']))
    
    if(!is.na(milestones[i, 'apv'])) {
      milestones[i, 'hpv'] <- growthfd.evaluate(milestones[i, 'apv'], fit$par, model)
    }

    if(!is.na(milestones[i, 'atf'])) {
      milestones[i, 'htf'] <- growthfd.evaluate(milestones[i, 'atf'], fit$par, model)
    }
    
    # Evaluation of fits and residuals
    msk.na <- id.na == ids[i]
    
    f <- growthfd.evaluate(x.na[msk.na], fit$par, model)
    r <- f - y.na[msk.na]
    fitted <- rbind(fitted, data.frame('id'=id.na[msk.na], 'fitted'=f, 'residuals'=r))
    
    # Evaluations of stature, velocity and acceleration
    stature[i,] <- growthfd.evaluate(sampling, fit$par, model)
    velocity[i,] <- growthfd.evaluate(sampling, fit$par, model, 1)
    acceleration[i,] <- growthfd.evaluate(sampling, fit$par, model, 2)
  }
  
  colnames(fitted) <- c('id', 'fitted', 'residuum')
  return(list('ids' = ids, 'scores' = scores, 'milestones' = milestones, 'fitted' = fitted, 'stature' = stature, 'velocity' = velocity, 'acceleration' = acceleration, 'sampling' = sampling))
}


#' Compute apv of model instance
#'
#' This function computes apv related to the certain instance of the model
#' described by the given parameters.
#'
#' @param model FPCA growth model
#' @param par Params of the model, corresponding to some individual
#' @return Age of maximum growth velocity
#' @export
growthfd.apv <- function(model, par) {
  x <- seq(7, 18, 0.05)
  y <- growthfd.evaluate(x, par, model, deriv=1)
  xyi <- data.frame(x, y)
  
  peak.xyi <- sitar::getPeak(xyi)
  return(peak.xyi[1])
}

#' Plot a velocity curve registered at apv
#'
#' This function plots a velocity curve, registered at population (model) apv
#' in comparison with the mean curve.
#'
#' @param model FPCA growth model
#' @param par Params of the model, corresponding to the individual
#' @example man/examples/plot.ApvRegVelocity.R
#' @return Velocity at apv plot
#' @export
growthfd.plot.ApvRegVelocity <- function(model, par) {
  meanPar <- rep(0, length(par))
  pApv <- growthfd.apv(model, meanPar)
  iApv <- growthfd.apv(model, par)
  
  message(sprintf('Population apv=%f, individual apv=%f', pApv, iApv))
  
  wbasisLM <- fda::create.bspline.basis(c(0,18), 4, 3, c(0, pApv, 18))
  WfdLM <- fda::fd(matrix(0,4,1), wbasisLM)
  WfdParLM <- fda::fdPar(WfdLM, 1, 1e-12)
  
  fd <- growthfd.std(par, model)
  reg <- fda::landmarkreg(fd, iApv, pApv, WfdParLM, FALSE)
  
  x <- seq(0,18,0.05)
  n <- length(x)
  yi <- fda::eval.fd(x, fda::register.newfd(fd, reg$warpfd), 1)
  ym <- growthfd.evaluate(x, meanPar, model, deriv=1)

  data <- data.frame(id=rep('Mean', n), x=x-pApv, y=ym)
  data <- rbind(data, data.frame(id=rep('Individual', n), x=x-pApv, y=yi))
  
  p <- ggplot2::ggplot(data=data, ggplot2::aes(x=x, y=mean, group=id, color=id)) +
    ggplot2::geom_line(ggplot2::aes(linetype=id)) +
    ggplot2::xlim(-3, 3) + 
    ggplot2::ylim(0, 13) +
    ggplot2::xlab('Years before and after time of maximum velocity') +
    ggplot2::ylab('Height velocity (cm/yr)') +
    ggplot2::scale_color_brewer(palette="Paired") +
    ggplot2::theme_minimal()
  
  return(p)
}

