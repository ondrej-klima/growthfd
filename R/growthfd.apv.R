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
