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
growthfd.bgs.gather <- function(data, prefix, age=NULL) {
  if(is.null(age)) {
    age <- growthfd.bgs.measurementsAge()
  }
  n <- sprintf('%s%.2d%.2d', rep(prefix, length(age)), floor(age), (age%%1)*1e2)
  s <- subset(data, select=c('id', 'sex', n)) %>% gather(agecat, value, n)
  s <- s[order(s$id),]
  s$age <- rep(age, dim(data)[1])
  return(s)
}