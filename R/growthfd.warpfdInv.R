#' Inverse time warping function
#' 
#' This function returns the *inverse* time warping function corresponding 
#' to supplied model and particular parameters.
#' 
#' @param par Parameters of the model
#' @param model FPCA growth model  
#' @import fda
#' @export
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
