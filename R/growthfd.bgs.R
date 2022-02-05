#' List ids of individuals to be dropped from height modeling
#'
#' This function returns a vector containing ids of individuals with incomplete 
#' stature measurements.
#'
#' @return Vector of ids
#' @export
growthfd.bgs.dropoutsIds.Height <- function() {
  return(c(
    0,2,3,4,
    11,12,14,18,
    25,28,
    44,47,49,
    58,
    70,71,78,
    83,85,86,87,88,
    91,93,96,
    100,104,105,
    116,
    120,123,124,127,
    130,133,
    144,148,149,
    152,155,156,157,
    166,169,
    174,
    182,183,184,186,
    191,194,
    202,207,
    213,214,218,219,
    222,228,
    231,234,236,
    240,243,244,249,
    251,253,
    261,265,
    275,278,
    280,283,285,
    290,292,294,297,
    310,311,318,319,
    322,324,
    332,
    344,347,
    355,356,
    365,
    372,
    388,
    394,396,397,
    401,402,
    415,418,
    420,425,
    432,
    440,443,
    457,
    461,462,464,
    471,
    487,
    516,
    521,522,523,
    531,538
    )
  )
}


#' List ages of measurements
#' 
#' This function returns a vector of ages when the measurements were performed.
#' 
#' @return Vector of ages
#' @export 
growthfd.bgs.measurementsAge <- function() {
  return(c(0, 0.25, 0.5, 0.75, seq(1, 18, 0.5)))
}

#' Gather selected columns
#' 
#' Selects columns with given prefix for supplied ages and gathers them into 
#' a matrix together with id and sex data
#' 
#' @param data BGS data
#' @param prefix Columns prefix
#' @param age Vector containing ages (optional)
#' @return Gathered data
#' @importFrom dplyr %>%
#' @export
growthfd.bgs.gather <- function(data, prefix='vysk', age=NULL) {
  if(is.null(age)) {
    age <- growthfd.bgs.measurementsAge()
  }
  n <- sprintf('%s%.2d%.2d', rep(prefix, length(age)), floor(age), (age%%1)*1e2)
  s <- subset(data, select=c('id', 'sex', n)) %>% tidyr::gather(agecat, value, n)
  s <- s[order(s$id),]
  s$age <- rep(age, dim(data)[1])
  return(s)
}

#' Estimate NA values
#' 
#' Interpolates missing values using spline method. Interpolates data 
#' from the 'value' column and join them as the 'valuei' column.
#'
#' @param gatheredData Data in gathered form
#' @return Interpolated data
#' @export
growthfd.bgs.interpolateNAs <- function(gatheredData) {
  ids <- unique(gatheredData$id)
  gatheredData$valuei <- rep(NA, length(gatheredData$value))
  for(id in ids) {
    rows <- gatheredData$id == id
    value <- gatheredData$value[rows]
    gatheredData$valuei[rows] <- imputeTS::na_interpolation(value, option = "stine")
  }
  return(gatheredData)
}

#' Resample the data 
#' 
#' Resamples the data without NA values to fine grid.
#' 
#' @param interpolatedData Data to be resampled.
#' @return Resampled data
growthfd.bgs.resample <- function(interpolatedData) {
  result <- list()
  ids <- unique(interpolatedData$id)
  for(id in ids) {
    rows <- interpolatedData$id == id
    br <- zoo::zoo(interpolatedData$valuei[rows], interpolatedData$age[rows])
    t <- time(br)
    t.4iy <- seq(from = min(t), to=max(t),by=0.05)
    dummy <- zoo::zoo(,t.4iy)
    br.interpolated <- merge(br,dummy,all=TRUE)
    br.spl <- zoo::na.spline(br.interpolated, method = "natural")
    w.spl <- time(br.spl)
    Tset.spl <- t.4iy 
    brs.spl <- br.spl[w.spl %in% Tset.spl]
    age <- as.vector(time(brs.spl))
    result[[id]] <- cbind('id'=rep(id,length(age)),'age'=age,'value'=as.vector(brs.spl))
  }
  return(do.call(rbind, result))
}

#' Fit the monotone spline
#'
#' This function fit the monotone splines to the data.
#' 
#' @param resampledData Data to be interpolated by monotone fda splines
#' @return Object with fitted splines
#' @example man/examples/train.R
#' @export
growthfd.bgs.smooth <- function(resampledData, monotone=T, norder=6, Lfdobj=3, lambda=5e-2) {
  age <- unique(resampled[,'age'])    
  values <- resampledData[,'value']
  ncases <- length(unique(resampledData[,'id']))
  dim(values) <- c(length(age), ncases)
  
  rng <- c(min(age),max(age))
  
  nage <- length(age)
  nbasis <- nage + norder - 2
  wbasis <- fda::create.bspline.basis(rng, nbasis, norder, age)
  cvec0 <- matrix(0,nbasis,ncases)
  Wfd0 <- fda::fd(cvec0, wbasis)
  growfdPar <- fda::fdPar(Wfd0, Lfdobj, lambda)
  wgt <- rep(1, nage)
  
  return(if(monotone) {
    fda::smooth.monotone(age, values, growfdPar, wgt, conv=0.1)
  } else {
    fda::smooth.basis(age, values, growfdPar, wgt)
  })
}

#' Evaluate monotone fda splines
#'
#' This function evaluates the set of monotone splines.
#'
#' @param fda Fda object
#' @param age Vector of ages to be evaluated
#' @param deriv Derivative of the curve
#' @return Matrix of evaluated points
#' @export
growthfd.bgs.evalMonotone <- function(fda, age, deriv=0) {
  n <- length(age)
  m <- dim(fda$beta)[2]
  result <- matrix(NA, n, m)
  for(i in seq(m)) {
    if(deriv == 0) {
      result[,i] <- fda$beta[1,i] + fda$beta[2,i]*fda::eval.monfd(age, fda$Wfdobj[i])
    } else {
      result[,i] <- fda$beta[2,i]*fda::eval.monfd(age, fda$Wfdobj[i], deriv)
    }
  }
  return(result)
}

#' Evaluate general fda splines
#'
#' This function evaluates non monotone fda splines.
#' 
#' @param fda Fda object
#' @param age Vector of ages to be evaluated
#' @param deriv Derivative of the curve
#' @return Matrix of evaluated points
#' @export
growthfd.bgs.eval <- function(fda, age, deriv=0) {
  return(fda::eval.fd(age, fda, deriv))
}

#' Find apvs on growth curves
#' 
#' This function finds ages of maximal growth velocity on the velocity curves.
#' 
#' @param age Vector of ages
#' @param velocity Matrix of velocity curves
#' @param ids Vector of individuals' ids
#' @param limits List of limits
#' @return Vector of apv values
#' @export
growthfd.bgs.apvs <- function(age, velocity) {
  n <- dim(velocity)[2]
  result <- rep(NA, n)
  for(col in seq(n)) {
    xyi <- data.frame(x=age, y=velocity[,col])
    result[col] <- sitar::getPeak(xyi)[1]
  }
  return(result)
}

#' Plot all individual curves to pdf
#' 
#' Plots value, velocity and acceleration curves together with apvs
#' and measured data points into separate figure for each individual. All plots
#' are stored into a single pdf file, one figure per page.
#' 
#' @param age Vector of ages
#' @param ids Vector containing ids
#' @param apvs Vector containing apv for each individual
#' @param values Matrix with value curves
#' @param vel Matrix with velocity curves
#' @param values Matrix with acceleration curves
#' @param data Matrix with original data points
#' @param filename File name of the output pdf
#' @export
growthfd.bgs.plotIndividuals <- function(age, ids, apvs, values, vel, acc, data, filename="plots.pdf") {
  n <- length(ids)
  pdf(filename)
  for(i in seq(n)) {
    df <- data.frame('age'=age, 'v'=values[,i], 'vel'=vel[,i], 'acc'=acc[,i])
    rows <- data$id == ids[i]
    dfPts <- data.frame('age' = data$age[rows], 'v' = data$value[rows])
    p <- ggplot2::ggplot(data = df) +
      ggplot2::geom_line(ggplot2::aes(x=age,y=v)) +
      ggplot2::geom_line(ggplot2::aes(x=age,y=acc)) +      
      ggplot2::geom_line(ggplot2::aes(x=age,y=vel)) +
      ggplot2::geom_vline(xintercept = apvs[i]) + 
      ggplot2::geom_point(data=dfPts, ggplot2::aes(x=age, y=v)) +
      ggplot2::ggtitle(ids[i])
    print(p)
  }
  dev.off()
}

#' Plot curves in one figure
#' 
#' Plots all curves from given matrix into a single figure.
#' 
#' @param age Vector of ages
#' @param values Matrix containing curves as columns
#' @param xlimit Limits for the x axis
#' @param ylimit Limits for the y axis
#' @return GGPlot2 plot
#' @export
growthfd.bgs.plotAll <- function(age, values, xlimit=NULL, ylimit=NULL) {
  d <- dim(values)
  dim(values) <- c(prod(dim(values)), 1)
  df <- data.frame(age=rep(age, d[2]), values=values, i=factor(rep(seq(d[2]), each=d[1])))
  result <- ggplot2::ggplot(data = df) + 
    ggplot2::geom_line(ggplot2::aes(x=age,y=values,group=i,colour=i),show.legend = FALSE)
  if(!is.null(xlimit)) {
    result <- result + ggplot2::xlim(xlimit)
  }
  if(!is.null(ylimit)) {
    result <- result + ggplot2::ylim(ylimit)
  }
  return(result)
}

#' Register curves to the apvs
#' 
#' Calculates the time warping functions with respect to the supplied apvs 
#' and the growth curves for the final refinement.
#'
#' @param fdaObject Fda object contaning the curves
#' @param apvs Vector containing the respective apvs
#' @return FDA object containing the time warping functions
#' @export
growthfd.bgs.registerCurvesToApvs <- function(fdaObject, apvs) {
  hgtfhatfd = fdaObject$yhatfd;
  accelfdUN = fda::deriv.fd(hgtfhatfd, 2)
  accelmeanfdUN = fda::mean.fd(accelfdUN)
  PGSctr=apvs
  PGSctrmean = mean(PGSctr)
  ncases <- length(apvs)
  
  wbasisLM = fda::create.bspline.basis(c(0,18), 4, 3, c(0,PGSctrmean,18))
  WfdLM = fda::fd(matrix(0,4,1),wbasisLM)
  WfdParLM = fda::fdPar(WfdLM,1,1e-12)
  
  regListLM = fda::landmarkreg(accelfdUN, PGSctr, PGSctrmean, WfdParLM, TRUE)
  
  accelfdLM = regListLM$regfd
  accelmeanfdLM = fda::mean.fd(accelfdLM)
  warpfdLM = regListLM$warpfd
  WfdLM = regListLM$Wfd
  
  wbasisCR = fda::create.bspline.basis(c(0,18), 15, 5)
  Wfd0CR = fda::fd(matrix(0,15,ncases),wbasisCR)
  WfdParCR = fda::fdPar(Wfd0CR, 1, 1)
  regList = fda::register.fd(accelmeanfdLM,accelfdLM, WfdParCR)
  warpfdCR = regList$warpfd
  
  return(fda::register.newfd(warpfdLM, warpfdCR))
}

#' Compute inverse time-warping functions
#' 
#' Computes inverse for given functions.
#' 
#' @param age Vector of ages
#' @param tw Fda object containing the time-warping functions
#' @return Fda object containing the inverse functions
#' @import splines 
#' @export 
growthfd.bgs.invertTw <- function(age, tw) {
  values=fda::eval.fd(age, tw)
  
  values[1,] = 0
  values[361,] = 18
  rng <-c(0,18)
  
  nage <- dim(values)[1]
  ncases <- dim(values)[2]
  norder <- 6
  nbasis <- nage + norder - 2
  Lfdobj <- 3 
  lambda <- 5e-2
  
  d <- matrix(nrow = nage, ncol = ncases)
  for(i in 1:ncases) {
    wbasis <- fda::create.bspline.basis(rangeval=rng, nbasis=nbasis, norder=norder, breaks=values[,i])
    cvec0 <- matrix(0,nbasis,1)
    Wfd0 <- fda::fd(cvec0, wbasis)
    growfdPar <- fda::fdPar(Wfd0, Lfdobj, lambda)
    
    d[,i] <- fda::eval.fd(age, smooth.basis(values[,i], age, growfdPar)$fd)
  }

  wbasis <- fda::create.bspline.basis(rangeval=rng, nbasis=nbasis, norder=norder, breaks=age)
  Wfd0 <- fda::fd(cvec0, wbasis)
  growfdPar <- fda::fdPar(Wfd0, Lfdobj, lambda)
  return(fda::smooth.basis(age, d, growfdPar))
}

#' Create FPCA growth model
#'
#' Creates FPCA growth model from fda objects of growth functions registered
#' on apv and inverse time warping functions.
#'
#' @param ampitude Fda object of registered growth functions
#' @param itw Fda object of inverse time warping functions
#' @param nharm Number of harmonic functions for each fpca
#' @return FPCA growth model
#' @export
growthfd.bgs.model <- function(amplitude, itw, nharm=6) {
  model<-list();
  model$growthfpca <- fda::pca.fd(amplitude,nharm=nharm);
  model$warpfpca <- fda::pca.fd(itw, nharm=nharm);
  model$scores.elements <- c(
    dim(model$warpfpca$scores)[2], 
    dim(model$growthfpca$scores)[2]
  )
  return(model)
}
