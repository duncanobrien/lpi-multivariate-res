#' Multivariate Jacobian Index Estimated From Multivariate Autocorrelation Matrix
#'
#' Estimate the dominant Jacobian eigenvalue of a multivariate time series using autocorrelated stochastic differential equations
#'
#' @param data Numeric matrix with time in first column and species abundance in the second
#' @param winsize Numeric. Defines the window size of the rolling window as a percentage of the time series length.
#' @param scale Boolean. Should data be scaled prior to estimating the Jacobian.
#' @param dt Numeric An appropriate time step
#' 
#' @returns A dataframe where the first column is last time index of the window and the second column is the estimated index value. A value <1.0 indicates stability, a value >1.0 indicates instability.
#'
#' @examples
#' #Load the multivariate simulated
#' #dataset `simTransComms`
#'
#' data(simTransComms)
#'
#'#Subset the second community prior to the transition
#'
#' pre_simTransComms <- subset(simTransComms$community2,time < inflection_pt)
#'
#' #Estimate the univariate stability index for the first species in
#' #the second community
#'
#' egarJ <- multiAR(data = pre_simTransComms[,2:3],
#' winsize = 25, dt = 1)
#'
#' @export
#' @source Williamson and Lenton (2015). Detection of bifurcations in noisy coupled systems from multiple time series. Chaos, 25, 036407

multiAR <- function(data, scale = TRUE, dt = 1, method = c("rolling","expanding"), 
                    winsize = 50, burn_in = 7, tail.direction = "one.tailed", threshold = 2){
  
  data <- as.data.frame(data)
  
  if(method == "rolling"){
  window <- round(dim(data)[1] * winsize/100)
  
  out <- lapply(1:(dim(data)[1]-window+1), function(i){
    
    sub_data <- data[i:(i+window-1),]
    if(isTRUE(scale)){
    sub_data <- sapply(sub_data[,-1],FUN = function(x){return(c(scale(x)))} )
    }
    mAr_mod <- mAr::mAr.est(as.matrix(sub_data[,-1]), 1)
    yy <- eigen(mAr_mod$AHat)
    jac_eig <- (1/dt)*(log(abs(yy$values)))
    return(cbind("time" = data[(i+window-1),1],
                 "multiAR" = max(jac_eig)))
    })
  raw <- as.data.frame(do.call("rbind", out))
  
  out.cor <- data.frame("multiAR" = tryCatch({cor.test(as.numeric(raw$time), raw$multiAR, alternative = c("two.sided"), method = c("kendall"), conf.level = 0.95,na.action = na.omit)$estimate}, error = function(e){ warning("Correlation coefficents not returned as too few observations"); return(NA)}))
  out <- list("raw" = raw, "cor" = out.cor)
  
  }
  
  if(method == "expanding"){
    RES<-list()
    roll.ar <-NULL
    
    for(i in (burn_in):dim(data)[1]){
     
     sub_data <- data[1:i,]
     
     if(isTRUE(scale)){
       sub_data <- sapply(sub_data[,-1],FUN = function(x){return(c(scale(x)))} )
     }
     
     mAr_mod <- mAr::mAr.est(as.matrix(sub_data[,-1]), 1)
     yy <- eigen(mAr_mod$AHat)
     jac_eig <- (1/dt)*(log(abs(yy$values)))
     
     roll.ar[[i]] <- max(jac_eig)
     multi.ar<- (max(jac_eig)-mean(unlist(roll.ar), na.rm=TRUE))/sd(unlist(roll.ar), na.rm = TRUE)
     
     RES[[i]] <- data.frame("time" = data[i,1],
                            "meanAR" = multi.ar)
     }
    
    results <- do.call("rbind", RES)
    
    out<-data.frame(results) %>%
      tidyr::pivot_longer(-c("time"), names_to = "metric.code",values_to = "metric.score") %>%
      dplyr::group_by(.data$metric.code) %>% dplyr::arrange(.data$time,.by_group = TRUE) %>%
      dplyr::mutate(rolling.mean = rolling_mean(.data$metric.score),
                    rolling.sd = rolling_sd(.data$metric.score))
    out$threshold.crossed<-NA
    
    if(tail.direction == "two.tailed"){
      out$threshold.crossed[which(out$metric.score>(out$rolling.mean+(threshold*out$rolling.sd))|out$metric.score<(out$rolling.mean-(threshold*out$rolling.sd)))]<-1
      out$threshold.crossed[is.na(out$threshold.crossed)] <- 0
    }else{
      out$threshold.crossed[which(out$metric.score>(out$rolling.mean+(threshold*out$rolling.sd)))]<-1
      out$threshold.crossed[is.na(out$threshold.crossed)] <- 0
    }
    
    out$str<-(out$metric.score-out$rolling.mean)/out$rolling.sd
    }
  
  return(out)
  }

rolling_mean <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- mean(x[1:i], na.rm=T);
  }
  return(result);
}

rolling_sd <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- sd(x[1:i], na.rm=T);
  }
  return(result);
}
