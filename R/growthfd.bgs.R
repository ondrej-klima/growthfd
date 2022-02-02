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
#' @export
growthfd.bgs.gather <- function(data, prefix='vysk', age=NULL) {
  if(is.null(age)) {
    age <- growthfd.bgs.measurementsAge()
  }
  n <- sprintf('%s%.2d%.2d', rep(prefix, length(age)), floor(age), (age%%1)*1e2)
  s <- subset(data, select=c('id', 'sex', n)) %>% gather(agecat, value, n)
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
    gatheredData$valuei[rows] <- na_interpolation(value, option = "stine")
  }
  return(gatheredData)
}

#' Resample the data 
#' 
#' Resample the data without NA values to fine grid.
#' 
#' @param interpolatedData Data to be resampled.
#' @return Resampled data
growthfd.bgs.resample <- function(interpolatedData) {
  result <- list()
  ids <- unique(interpolatedData$id)
  for(id in ids) {
    rows <- interpolatedData$id == id
    br <- zoo(interpolatedData$valuei[rows], interpolatedData$age[rows])
    t <- time(br)
    t.4iy <- seq(from = min(t), to=max(t),by=0.05)
    dummy <- zoo(,t.4iy)
    br.interpolated <- merge(br,dummy,all=TRUE)
    br.spl <- na.spline(br.interpolated, method = "natural")
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
  wbasis <- create.bspline.basis(rng, nbasis, norder, age)
  cvec0 <- matrix(0,nbasis,ncases)
  Wfd0 <- fd(cvec0, wbasis)
  growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
  wgt <- rep(1, nage)
  
  return(if(monotone) {
    smooth.monotone(age, values, growfdPar, wgt, conv=0.1)
  } else {
    smooth.basis(age, values, growfdPar, wgt)
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
  return(if(deriv == 0) {
    fda$beta[1,] + fda$beta[2,]*eval.monfd(age, fda$Wfdobj)
  } else {
    fda$beta[2,]*eval.monfd(age, fda$Wfdobj, deriv)
  })
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
  return(eval.fd(age, fda, deriv))
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
growthfd.bgs.apvs <- function(age, velocity, ids, limits) {
  n <- dim(velocity)[2]
  result <- rep(NA, n)
  for(col in seq(n)) {
    xyi <- data.frame(x=age, y=velocity[,col])
    result[col] <- sitar::getPeak(xyi)[1]
  }
  return(result)
}

