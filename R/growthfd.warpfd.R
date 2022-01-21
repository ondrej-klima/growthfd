#' Time warping function
#' 
#' This function returns the time warping function corresponding 
#' to supplied model and particular parameters.
#' 
#' @param par Parameters of the model
#' @param model FPCA growth model  
#' @import fda
#' @export
growthfd.warpfd <- function(par, model) {
  nwarpscores <- model$scores.elements[1];
  warpscores <- par[1:nwarpscores];
  
  warpfd <- model$warpfpca$meanfd;
  for(i in 1:nwarpscores) {
    warpfd <- warpfd + warpscores[i]*sqrt(model$warpfpca$values[i])*model$warpfpca$harmonics[i];
  }
  
  return(warpfd)
}
