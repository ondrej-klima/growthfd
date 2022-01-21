#' Prepare data for velocity boxplots registered on apv
#' 
#' This function prepares data for boxplots in time of measurements, 
#' after registration of the individual on population apv.
#' 
#' @param model Model
#' @param par Parameters of the model fitted to the measurements
#' @param ages Ages of measurements points
#' @param rndn Count of random curves to be evaluated
#' @return Data frames for population and the individual
#' @export
growthfd.RegVelocities <- function(model, par, ages, rndn = 0, verbose=F) {
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
    invisible(
      capture.output(
        reg <- fda::landmarkreg(fd, iApv, pApv, WfdParLM, FALSE)
        )
      )
  }
    
  yi <- fda::eval.fd(ages, fda::register.newfd(fd, reg$warpfd), 1)
  
  if(rndn == 0) {
    growthScores <- growthfd.modelPars(model)[,7:12]
    n <- nrow(growthScores)
    scores <- cbind(matrix(0,n,6), growthScores);
    if(verbose) {
      message(sprintf('Using %d individuals from the model.', n))
    }
  } else {
    scores <- cbind(matrix(0,rndn,6), MASS::mvrnorm(rndn,rep(0,6),diag(6)));
    if(verbose) {
      message(sprintf('Using %d random individuals.', rndn))
    }
    n <- rndn;
  }
  
  y <- matrix(NA, n, length(ages))
  for(i in seq(n)) {
    if(verbose) {
      message(sprintf('Evaluating curve %d..', i))
    }
    y[i,] <- growthfd.evaluate(ages, scores[i,], model, deriv = 1)
  }
  
  return(list(
    'individual' = data.frame(
      x=as.factor(round(ages-pApv, digits=2)), 
      y = yi), 
    'population' = data.frame(
      x=as.factor(round(rep(ages, n)-pApv, digits=2)), 
      y=c(t(y)))
    )
  )
}
