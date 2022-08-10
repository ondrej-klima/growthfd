#' Fit a FPCA Growth Curve Model to a population
#'
#' This function fits a model to the given measured data of a population.
#'
#' @param model FPCA growth model to be fitted
#' @param data Data frame containing age, height and id of individuals
#' @param x Age at measured data points
#' @param y Height at measured data points
#' @param id Corresponding individual's id at measured data points
#' @param verbose Verbosity
#' @param bounds Limitation of the interval for milestones estimation, 'negative' or 'inverse'
#' @param filename File name for saving results after each individual
#' @param startFromId Start the evaluation from this id
#' @param parallel (Experimental) Parallel evaluation of the model fitting
#' @param scores.filename File name for continuous saving of the scores
#' @return List containing individuals id and model 
#' @example man/examples/main.R
#' @importFrom foreach %dopar%
#' @export
growthfd <- function(data, 
                     x, 
                     y, 
                     id, 
                     model, 
                     verbose=1, 
                     bounds='negative', 
                     filename = '', 
                     startFromId = NULL, 
                     parallel = F, 
                     scores.filename = 'parallel.txt') {
  
  mcall <- match.call()
  x.na <- as.numeric(eval(mcall$x, data))
  y.na <- as.numeric(eval(mcall$y, data))
  id.na <- as.factor(eval(mcall$id, data))
    
  msk <- !is.na(x.na) & !is.na(y.na) & x.na <= 18 & x.na >= 8
  x <- x.na[msk]
  y <- y.na[msk]
  id <- id.na[msk]
  
  msk <- !is.na(x.na) & x.na <= 18 & x.na >= 8
  x.na <- x.na[msk]
  y.na <- y.na[msk]
  id.na <- id.na[msk]

  #ids <- levels(id)
  ids <- unique(id)
  n <- length(ids)
  
  scores <- matrix(NA, n, sum(model$scores.elements))
  milestones <- matrix(NA, n, 6)
  colnames(milestones) <- c('apv', 'vpv', 'hpv', 'atf', 'vtf', 'htf')
  w.m <- matrix(NA, n, 2)
  colnames(w.m) <- c('apv', 'atf')
  fitted <- data.frame('id'=factor(), 'fitted'=double(), 'residuals'=double(),
                       'age' = double(), 'apv' = double())
  
  # sampling <- seq(0, 18, 0.25)
  sampling <- seq(8, 19.5, 0.25)
  m <- length(sampling)
  stature <- matrix(NA, n, m)
  velocity <- matrix(NA, n, m)
  acceleration <- matrix(NA, n, m)
  
  # Determine model apv and atf
  m_x <- seq(8, 18, 0.05)
  m_y <- growthfd.evaluate(m_x, 
                           rep(0, sum(model$scores.elements)), 
                           model, 
                           deriv=1)
  xyi <- data.frame('x'=m_x, 'y'=m_y)
  m_apv <- peak.xyi <- sitar::getPeak(xyi)[1]
  m_atf <- sitar::getTakeoff(xyi)[1]
  message(sprintf('Model apv=%f, atf=%f', m_apv, m_atf))
  
  fromId <- if(is.null(startFromId)) { 
    1 
  } 
  else { 
    match(startFromId, ids) 
  }
  
  start_t <- proc.time()
  if(!parallel) {
    for(i in seq(fromId, length(ids))) {   
      msk <- id == ids[i]
  
      if(!any(msk)) {
        #message('Skipping due to no data..')
        next
      }
      
      if(verbose > 0) {
        message(sprintf("Processing individual with id '%s' (%d/%d), containing %d measurements\n", ids[i], i, n, sum(msk)))
      }
      # message(x[msk])
      # message(y[msk])
      
      fit <- growthfd.fit(model, x[msk], y[msk], verbose)
      scores[i,] <- fit$par
      
      if(filename != '') {
        growthfd.result <- list('ids' = ids, 'scores' = scores, 'lastProcessedId' = ids[i])
        save(growthfd.result, file = filename)
      }
    }
  }
  else {
    if(scores.filename != '') {
      cat('i,id,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12\n', file=scores.filename)
    }
    lock<-tempfile()
    cores=parallel::detectCores()
    cl<-parallel::makeCluster(cores-1)
    doParallel::registerDoParallel(cl)
    scores <- foreach::foreach(i=seq(fromId, length(ids)), .combine=rbind, .packages=c('growthfd'), .verbose=T) %dopar% {
      msk <- id == ids[i]
      fit <- rep(NA, 12)
      if(any(msk)) {
        if(verbose > 0) {
          sink(file=sprintf('%d.txt', Sys.getpid()), append = TRUE, type = c("output", "message"))
          cat(sprintf("Processing individual with id '%s' (%d/%d), containing %d measurements\n", ids[i], i, n, sum(msk)))
        }
        fit <- growthfd.fit(model, x[msk], y[msk], verbose)$par
        if(verbose > 0) {
          if(scores.filename != '') {
            cat('Writing resulting scores to file...\n')
          }
          sink()
        }
        if(scores.filename != '') {
          locked <- flock::lock(lock)
          write.table(t(c(i, ids[i], fit)), file=scores.filename, col.names = F, row.names = F, sep=",", append = TRUE)
          flock::unlock(locked)
        }
      }
      fit
    }
    parallel::stopCluster(cl)
    if(filename != '') {
      growthfd.result <- list('ids' = ids, 'scores' = scores, 'lastProcessedId' = ids[length(ids)])
      save(growthfd.result, file = filename)
    }
    closeAllConnections()
  }
  fit_t <- proc.time() - start_t
  message('Fitting time:\n')
  print(fit_t)
  
  for(i in seq(fromId, length(ids))) { 
    if(!any(id == ids[i])) {
      #message('Skipping due to no data..')
      next
    }
    
    if(bounds == 'negative') {
      warpfd <- growthfd.warpfd(-scores[i,], model)
    }
    else {
      warpfd <- growthfd.warpfdInv(scores[i,], model)
    }
    w_m <- fda::eval.fd(c(m_apv, m_atf), warpfd)
    #w_apv <- 2*m_apv - w_m[1]
    #w_atf <- 2*m_atf - w_m[2]
    w_apv <- w_m[1]
    w_atf <- w_m[2]
    message(sprintf('Warped apv=%f, atf=%f', w_apv, w_atf))
    w.m[i, 'apv'] <- w_apv
    w.m[i, 'atf'] <- w_atf
    
    # Computation of growth milestones
    bound <- 2
    # Upper limit
    upper <- w_apv + bound
    if(upper > 18) {
      upper <- 18
    }
    
    # apv
    m_x_apv <- seq(w_apv-bound, upper, 0.05)
    # m_x_apv <- seq(7, 18, 0.05)
    m_y_apv <- growthfd.evaluate(m_x_apv, scores[i,], model, deriv=1)
    xyi_apv <- data.frame('x'=m_x_apv, 'y'=m_y_apv)
    
    peak.xyi_apv <- sitar::getPeak(xyi_apv)
    milestones[i, 'apv'] <- peak.xyi_apv[1]
    milestones[i, 'vpv'] <- peak.xyi_apv[2]
    
    # atf
    #m_x_atf <- seq(w_atf-bound, upper, 0.05)
    m_x_atf <- seq(8, upper, 0.05)
    # m_x_atf <- seq(7, 18, 0.05)
    m_y_atf <- growthfd.evaluate(m_x_atf, scores[i,], model, deriv=1)
    xyi <- data.frame('x'=m_x_atf, 'y'=m_y_atf)
    
    
    takeoff.xyi <- sitar::getTakeoff(xyi)
    milestones[i, 'atf'] <- takeoff.xyi[1]
    milestones[i, 'vtf'] <- takeoff.xyi[2]
    
    message(sprintf('Refined apv=%f, atf=%f', milestones[i, 'apv'], milestones[i, 'atf']))
    
    if(!is.na(milestones[i, 'apv'])) {
      milestones[i, 'hpv'] <- growthfd.evaluate(milestones[i, 'apv'], scores[i,], model)
    }

    if(!is.na(milestones[i, 'atf'])) {
      milestones[i, 'htf'] <- growthfd.evaluate(milestones[i, 'atf'], scores[i,], model)
    }
    
    
    # Evaluation of fits and residuals
    msk.na <- id.na == ids[i]
    
    f <- growthfd.evaluate(x.na[msk.na], scores[i,], model)
    r <- f - y.na[msk.na]
    fitted <- rbind(fitted, data.frame('id'=id.na[msk.na], 
                                       'fitted'=f, 
                                       'residuals'=r, 
                                       'age' = x.na[msk.na],
                                       'apv' = rep(milestones[i, 'apv'], length(f))))
    
    # Evaluations of stature, velocity and acceleration
    stature[i,] <- growthfd.evaluate(sampling, scores[i,], model)
    velocity[i,] <- growthfd.evaluate(sampling, scores[i,], model, 1)
    acceleration[i,] <- growthfd.evaluate(sampling, scores[i,], model, 2)
  }
  
  total_t <- proc.time() - start_t
  message('Total time:\n')
  print(total_t)
  
  
  colnames(fitted) <- c('id', 'fitted', 'residuum', 'age', 'apv')
  growthfd.result <-list('ids' = ids, 'scores' = scores, 'milestones' = milestones, 
       'fitted' = fitted, 'stature' = stature, 'velocity' = velocity, 
       'acceleration' = acceleration, 'sampling' = sampling, 'wm' = w.m,
       'lastProcessedId' = ids[i], 'time.fit' = fit_t, 'time.total' = total_t)
  
  if(filename != '') {
    save(growthfd.result, file = filename)
  }
  
  return(growthfd.result)
}
