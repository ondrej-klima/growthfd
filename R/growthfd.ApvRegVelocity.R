#' Register a velocity curve at population apv
#' 
#' This function registers a curve corresponding to the supplied parameters 
#' onto the population apv.
#' 
#' @param model FPCA growth model
#' @param par Params of the model, corresponding to the individual
#' @example man/examples/ApvRegVelocity.R
#' @return Velocity at apv data frame
#' @export
#'   
growthfd.ApvRegVelocity <- function(model, par, verbose = F) {
  meanPar <- rep(0, length(par))
  pApv <- growthfd.apv(model, meanPar)
  iApv <- growthfd.apv(model, par)
  
  if(verbose) {
    message(sprintf('Population apv=%f, individual apv=%f', pApv, iApv))
  }
  
  wbasisLM <- fda::create.bspline.basis(c(0,18), 4, 3, c(0, pApv, 18))
  WfdLM <- fda::fd(matrix(0,4,1), wbasisLM)
  WfdParLM <- fda::fdPar(WfdLM, 1, 1e-12)
  
  fd <- growthfd.std(par, model)
  
  if(verbose) {
    reg <- fda::landmarkreg(fd, iApv, pApv, WfdParLM, FALSE)
  } else {
    invisible(capture.output(reg <- fda::landmarkreg(fd, iApv, pApv, WfdParLM, FALSE)))
  }
  
  x <- seq(0,18,0.05)
  n <- length(x)
  yi <- fda::eval.fd(x, fda::register.newfd(fd, reg$warpfd), 1)
  ym <- growthfd.evaluate(x, meanPar, model, deriv=1)
  
  data <- data.frame(id=rep('Mean', n), x=x-pApv, y=ym)
  data <- rbind(data, data.frame(id=rep('Individual', n), x=x-pApv, y=yi))
}
