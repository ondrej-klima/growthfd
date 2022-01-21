#' Plot a velocity curve registered at apv
#'
#' This function plots a velocity curve, registered at population (model) apv
#' in comparison with the mean curve.
#'
#' @param model Data obtained using growthfd.ApvRegVelocity
#' @example man/examples/plot.ApvRegVelocity.R
#' @return Velocity at apv plot
#' @export
growthfd.plot.ApvRegVelocity <- function(data) {
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
