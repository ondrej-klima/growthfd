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
