#' Plot velocity boxplots registered on apv
#' 
#' This function plots boxplots in time of measurements, after registration 
#' of the individual on population apv.
#' 
#' @param populationData Data frame for population box plots
#' @param individualData Data frame for the individual
#' @return GGPlot2 plot
#' @example man/examples/plot.RegVelocities.R
#' @export
growthfd.plot.RegVelocities <- function(populationData, individualData) {
  p <- ggplot2::ggplot(data=populationData, ggplot2::aes(x=x, y=y)) + 
    ggplot2::geom_boxplot() +
    ggplot2::geom_line(data=individualData, ggplot2::aes(x=x, y=mean, group=1)) +
    ggplot2::geom_point(data=individualData, ggplot2::aes(x=x, y=mean, group=1), shape=8) +
    ggplot2::xlab('Years before and after time of maximum velocity') +
    ggplot2::ylab('Height velocity (cm/yr)') +
    ggplot2::scale_color_brewer(palette="Paired") + 
    ggplot2::theme_minimal() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p)
}
