#' Plot velocity boxplots registered on apv
#' 
#' This function plots boxplots in time of measurements, after registration 
#' of the individual on population apv.
#' 
#' @param model Model
#' @param par Parameters of the model fitted to the measurements
#' @param ages Ages of measurements points
#' @param rndn Count of random curves to be evaluated
#' @return GGPlot2 plot
#' @example man/examples/plot.RegVelocities.R
#' @export
growthfd.plot.RegVelocities <- function(model, par, ages, rndn = 0) {
  meanPar <- rep(0, length(par))
  pApv <- growthfd.apv(model, meanPar)
  iApv <- growthfd.apv(model, par)
  
  message(sprintf('Population apv=%f, individual apv=%f', pApv, iApv))
  
  wbasisLM <- fda::create.bspline.basis(c(0,18), 4, 3, c(0, pApv, 18))
  WfdLM <- fda::fd(matrix(0,4,1), wbasisLM)
  WfdParLM <- fda::fdPar(WfdLM, 1, 1e-12)
  
  fd <- growthfd.std(par, model)
  reg <- fda::landmarkreg(fd, iApv, pApv, WfdParLM, FALSE)
  
  yi <- fda::eval.fd(ages, fda::register.newfd(fd, reg$warpfd), 1)
  
  if(rndn == 0) {
    growthScores <- growthfd.modelPars(model)[,7:12]
    n <- nrow(growthScores)
    scores <- cbind(matrix(0,n,6), growthScores);
    message(sprintf('Using %d individuals from the model.', n))
  } else {
    scores <- cbind(matrix(0,rndn,6), MASS::mvrnorm(rndn,rep(0,6),diag(6)));
    message(sprintf('Using %d random individuals.', rndn))
    n <- rndn;
  }
  
  y <- matrix(NA, n, length(ages))
  for(i in seq(n)) {
    message(sprintf('Evaluating curve %d..', i))
    y[i,] <- growthfd.evaluate(ages, scores[i,], model, deriv = 1)
  }
  
  data <- data.frame(x=as.factor(round(rep(ages, n)-pApv, digits=2)), y=c(t(y)))
  datai <- data.frame(x=as.factor(round(ages-pApv, digits=2)), y = yi)
  
  p <- ggplot2::ggplot(data=data, ggplot2::aes(x=x, y=y)) + 
    ggplot2::geom_boxplot() +
    ggplot2::geom_line(data=datai, ggplot2::aes(x=x, y=mean, group=1)) +
    ggplot2::geom_point(data=datai, ggplot2::aes(x=x, y=mean, group=1), shape=8) +
    ggplot2::xlab('Years before and after time of maximum velocity') +
    ggplot2::ylab('Height velocity (cm/yr)') +
    ggplot2::scale_color_brewer(palette="Paired") + 
    ggplot2::theme_minimal() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p)
}
