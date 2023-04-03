source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/uni_smap_jacobian.R")
source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/uni_smapJI.R")

source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/mvi.R")

source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/GFisher.R")
source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/NFisherpdf.R")
source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/roundTO.R")

source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/smapJ_index.R")
source("/Users/ul20791/Desktop/Academia/PhD/Repositories/EWSmethods/EWSmethods/R/smap_jacobian.R")


parallel_multiJI <- function(dt,var,n_cores = 4,...){
  
  require(foreach)
  require(doFuture)
  require(progressr)
  
  doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  
  dt2 <- split(dt,by = var)
  pb <- progressr::progressor(steps = length(dt2))
  
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.export = c("multi_smap_jacobian","smap_jacobian_est","multiJI")) %dopar% {
    
    #out <- EWSmethods::smapJI(data = x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE],...)
    out <- multiJI(data = as.data.frame(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE]),...)
    out <- data.frame("var" = x[,var,with=F][1],out)
    names(out) <- c(var,"time","multiJI")
    
    pb()
    
    return(out)
  }
  
  return(data.table::rbindlist(dt2))
}

parallel_uniJI <- function(dt,var,n_cores = 4,...){
  
  require(foreach)
  require(doFuture)
  require(progressr)
  require(rEDM)
  doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  
  dt2 <- split(dt,by = var)
  pb <- progressr::progressor(steps = length(dt2))
  
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods","rEDM"),.export = c("uni_smap_jacobian","uni_smap_jacobian_est","uniJI")) %dopar% {
    
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
    tmp <- lapply(names(x)[grepl("spp_",names(x))],function(i){

      if(i == "spp_1"){
      #EWSmethods::uniJI(x[,c("time",i),with =FALSE],E = 5,winsize=50)
        uniJI(data= as.data.frame(x[,c("time",i),with =FALSE]),...)

        }else{
      #EWSmethods::uniJI(x[,c("time",i),with =FALSE],E = 5,winsize=50)[,-1]
          uniJI(as.data.frame(x[,c("time",i),with =FALSE]),...)[,-1]
          }
      })

    tmp <- do.call("cbind", tmp)
    
    out <- data.frame("var" = x[,var,with=F][1],
                      "time" = tmp[,1],
                      "mean_uniJI" = rowMeans(tmp[,-1],na.rm = TRUE),
                      "max_uniJI" = apply(tmp[,-1], 1, max,na.rm = TRUE))
    
    names(out) <- c(var,"time","mean_uniJI","max_uniJI")
    
    pb()
    
    return(out)
  }
  
  return(data.table::rbindlist(dt2))
}

parallel_FI <- function(dt,var,n_cores = 4,...){
  
  require(foreach)
  require(doFuture)
  require(progressr)
  
  doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  
  dt2 <- split(dt,by = var)
  pb <- progressr::progressor(steps = length(dt2))
  
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.export = c("FI","NFisherpdf","roundTO")) %dopar% {
    
    sost <- t(apply(x[,names(x)[grepl("spp_",names(x))],with =FALSE], MARGIN = 2, FUN = sd))
    #size of states. Transpose required to ensure a 1 x n matrix
    
    # out <- EWSmethods::FI(x[,"time",with =FALSE][[1]],as.data.frame(x[,names(x)[grepl("spp_",names(x))],with =FALSE]),sost = sost,winsize=50)
    # 
    # out <- data.frame("var" = x[,var,with=F][1],
    #              "time" = apply(out$t_win, MARGIN = 2, FUN = max),
    #              "FI" =out$FI)
    # names(out) <- c(var,"time","FI")
    # 
    
    #out <- EWSmethods::FI(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE],sost = sost,winsize=50)
    out <- FI(as.data.frame(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE]),sost = sost,...)
    out <- data.frame("var" = x[,var,with=F][1],out$FI)
    names(out) <- c(var,"time","FI")
    
    pb()
    
    return(out)
  }
  
  return(data.table::rbindlist(dt2))
}

parallel_mvi <- function(dt,var,n_cores = 4,...){
  
  require(foreach)
  require(doFuture)
  require(progressr)
  
  doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  
  dt2 <- split(dt,by = var)
  pb <- progressr::progressor(steps = length(dt2))
  
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","EWSmethods"),.export = c("mvi")) %dopar% {
    
    #out <- EWSmethods::mvi(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE],winsize=50)
    out <- mvi(as.data.frame(x[,c("time",names(x)[grepl("spp_",names(x))]),with =FALSE]),...)
    out <- data.frame("var" = x[,var,with=F][1],out)
    names(out) <- c(var,"time","mvi")
    
    pb()
    
    return(out)
  }
  
  return(data.table::rbindlist(dt2))
}
