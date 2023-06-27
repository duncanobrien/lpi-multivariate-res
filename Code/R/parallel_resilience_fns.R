
parallel_multiJI <- function(dt,var,n_cores = 4,sample_spp = TRUE,winsize = 50,...){
  
  require(foreach)
  require(doRNG)
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  dt2 <- split(dt,by = var)

  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.errorhandling = "remove") %dorng% {
      
    if(NROW(x) < 10){ #10 is the minimum viable time series
      out <- data.frame("var" = x[,var,with=F][1], seq(10 * winsize/100,10),NA)
      names(out) <- c(var,"time","multiJI")
    }else{
      window <- round(dim(x)[1] * winsize/100)
      spp_col <- names(x)[grepl("spp_",names(x))]
      
      #sum(x[[j]] == as.numeric(names(sort(table(x[[j]]),decreasing = T)[1])))< (window-1)
      
      ts_check <- sapply(spp_col,FUN = function(j){length(unique(x[[j]])) == 1 |  (sum(as.numeric(base::table(x[[j]]) %in% c(1,2))) == 1 & length(unique(x[[j]])) == 2) | any(rle(x[[j]])$lengths >= (window-3))})
      spp_col <- spp_col[!ts_check]
      n_spp_ts <- length(spp_col)
      
      if((window - 1) < n_spp_ts){
        #target_spp <- names(x[,spp_col,with = F])[colSums(x[,spp_col,with = F]) >0]
        
        target_spp <- spp_col[sapply(spp_col,function(spp){
          sum(x[,spp,with=F] == as.numeric(names(sort(table(x[,spp,with=F]),decreasing = T)[1])))< (window-1) })]
        
        if(isFALSE(sample_spp)){
          stop(paste0("Length of windowed time series (",window,") < number of nodes (" ,n_spp_ts,"). Increase winsize or allow sampling of species"))
        }
        if(isTRUE(sample_spp)){
          warning(paste0("Length of windowed time series (",window,") < number of nodes (" ,n_spp_ts,"). Species have been sampled to allow index calculation"))
          
          sampl_spp <- sample(target_spp,round(window-3),replace = F)
          # out <- tryCatch(EWSmethods::multiJI(data = na.omit(as.data.frame(x[,c("time",sampl_spp),with =FALSE])),winsize = winsize,...),
          #                 error = function(err){return(cbind(time = seq(dim(x)[1] * winsize/100,dim(x)[1]),multiJI = NA))})
          out <- tryCatch(EWSmethods::multiJI(data = na.omit(as.data.frame(x[,c("time",sampl_spp),with =FALSE])),winsize = winsize, scale = T),
                                           error = function(err){return(cbind(time = seq(dim(x)[1] * winsize/100,dim(x)[1]),multiJI = NA))})
        }
      # }else if(n_spp_ts < 3){
      #   out <- data.frame(seq(10 * winsize/100,10),NA)
        }else if(n_spp_ts < 2){
          out <- data.frame(seq(10 * winsize/100,10),NA)
      }else{
        out <- tryCatch(EWSmethods::multiJI(data = na.omit(as.data.frame(x[,c("time",spp_col),with =FALSE])),winsize = winsize, scale = T),
                        error = function(err){return(cbind(time = seq(dim(x)[1] * winsize/100,dim(x)[1]),multiJI = NA))})
        
      }
      out <- data.frame("var" = x[,var,with=F][1],out)
      names(out) <- c(var,"time","multiJI")
      
      #pb()
    }
    return(out)
  }
  print("\u2713 multiJI")
  parallel::stopCluster(cl)
  return(data.table::rbindlist(dt2))
}

parallel_uniJI <- function(dt,var,n_cores = 4,winsize = winsize, tau = 1 , E = 1){
  
  require(foreach)

  # doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  # future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  # 
  # progressr::handlers(global = TRUE)
  # progressr::handlers("cli")
  # 
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  dt2 <- split(dt,by = var)
  #pb <- progressr::progressor(steps = length(dt2))
  
  #dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods","rEDM"),.export = c("uni_smap_jacobian","uni_smap_jacobian_est","uniJI"),.errorhandling = "remove") %dopar% {
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods","rEDM"),.errorhandling = "remove") %dopar% {
      
    #nams <- names(x)[grepl("spp_",names(x))]
    # 
    # tmp <- matrix(NA, nrow = length(x[[1]])-round(length(x[[1]])*0.5)+1,ncol = length(nams)+1)
    # for(i in seq_along(nams)){
    #   
    #   if(i == 1){
    #     
    #     tt <-  uniJI(data= as.data.frame(x[,c("time",nams[i]),with =FALSE]),E = 5,tau = -10,winsize = 50)
    #     
    #     tmp[,1] <- tt[,1]
    #     tmp[,2] <- tt[,2]
    #   
    #   }else{
    #     tmp[,i+1] <- uniJI(data= as.data.frame(x[,c("time",nams[i]),with =FALSE]),E = 5,tau = -10,winsize = 50)[,-1]
    #   }
    #   
    # }
    tmp <- lapply(names(x)[grepl("spp_",names(x))],function(k){
      
      if(k == "spp_1"){
        #EWSmethods::uniJI(x[,c("time",i),with =FALSE],E = 5,winsize=50)
        EWSmethods::uniJI(data= as.data.frame(x[,c("time",k),with =FALSE]),winsize = winsize,E = E,tau = tau,scale = T)
        
      }else{
        #EWSmethods::uniJI(x[,c("time",i),with =FALSE],E = 5,winsize=50)[,-1]
        EWSmethods::uniJI(data = as.data.frame(x[,c("time",k),with =FALSE]),winsize = winsize,E = E,tau = tau,scale = T)[,-1]
      }
    })
    
    tmp <- do.call("cbind", tmp)
    #names(tmp) <- c("time",paste("uniJI",names(x)[grepl("spp_",names(x))],sep = "_"))
    
    out <- data.frame("var" = x[,var,with=F][1],
                      tmp,
                      "mean_uniJI" = rowMeans(tmp[,-1],na.rm = TRUE),
                      "max_uniJI" = apply(tmp[,-1], 1, max,na.rm = TRUE))
    
    names(out) <- c(var,"time",paste("uniJI",names(x)[grepl("spp_",names(x))],sep = "_"),"mean_uniJI","max_uniJI")
    
    #pb()
    
    return(out)
  }
  print("\u2713 uniJI")
  #closeAllConnections()
  parallel::stopCluster(cl)
  
  #future::plan(sequential)
  return(data.table::rbindlist(dt2))
}

parallel_FI <- function(dt,var,n_cores = 4,winsize = 50, TL = 90){
  
  require(foreach)

  # doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  # future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  # 
  # progressr::handlers(global = TRUE)
  # progressr::handlers("cli")
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  dt2 <- split(dt,by = var)
  #pb <- progressr::progressor(steps = length(dt2))
  
  #dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.export = c("FI","NFisherpdf","roundTO"),.errorhandling = 'remove') %dopar% {
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.errorhandling = 'remove') %dopar% {
      
    spp_col <- names(x)[grepl("spp_",names(x))]
    
    sampl_spp <- spp_col[!sapply(spp_col,function(spp){
      length(unique(x[[spp]])) == 1 |  (sum(as.numeric(base::table(x[[spp]]) == 1)) == 1 & length(unique(x[[spp]])) == 2)
    })]
    
    #sampl_spp <- sample(sampl_spp,5,replace = F)
    
    sost <- t(apply(x[,sampl_spp,with =FALSE], MARGIN = 2, FUN = sd))
    #size of states. Transpose required to ensure a 1 x n matrix
    
    # out <- EWSmethods::FI(x[,"time",with =FALSE][[1]],as.data.frame(x[,names(x)[grepl("spp_",names(x))],with =FALSE]),sost = sost,winsize=50)
    # 
    # out <- data.frame("var" = x[,var,with=F][1],
    #              "time" = apply(out$t_win, MARGIN = 2, FUN = max),
    #              "FI" =out$FI)
    # names(out) <- c(var,"time","FI")
    # 
    
    #out <- EWSmethods::FI(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE],sost = sost,winsize=50)
    out <- EWSmethods::FI(as.data.frame(x[,c("time",sampl_spp),with =FALSE]),sost = sost,winsize = winsize, TL = TL)
    out <- data.frame("var" = x[,var,with=F][1],out$FI)
    names(out) <- c(var,"time","FI")
    
    #pb(sprintf("x=%g"),x)
    
    return(out)
  }
  print("\u2713 FI")
  #closeAllConnections()
  parallel::stopCluster(cl)
  #future::plan(sequential)
  
  return(data.table::rbindlist(dt2))
}

parallel_mvi <- function(dt,var,n_cores = 4,winsize = 50,...){
  
  require(foreach)
  
  # doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  # future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  # 
  # progressr::handlers(global = TRUE)
  # progressr::handlers("cli")
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  dt2 <- split(dt,by = var)
  #pb <- progressr::progressor(steps = length(dt2))
  
  #dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.export = c("mvi"),.errorhandling = 'remove') %dopar% {
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.errorhandling = 'remove') %dopar% {
      
    #out <- EWSmethods::mvi(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE],winsize=50)
    out <- EWSmethods::mvi(as.data.frame(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE]),winsize = winsize)
    out <- data.frame("var" = x[,var,with=F][1],out)
    names(out) <- c(var,"time","mvi")
    
    #pb()
    
    return(out)
  }
  print("\u2713 mvi")
  #closeAllConnections()
  
  parallel::stopCluster(cl)
  #future::plan(sequential)
  return(data.table::rbindlist(dt2))
}
