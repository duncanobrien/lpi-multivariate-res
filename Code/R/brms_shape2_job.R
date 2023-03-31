brms_shape2 <- function(dt,var,n_cores = 4){
  
  require(foreach)
  require(brms)
  require(rstan)
  dt2 <- split(dt,by = var)

  cl <- parallel::makePSOCKcluster(floor(n_cores/2),outfile="Log.txt") # number of cores to use

  doSNOW::registerDoSNOW(cl)
  
  pb <- txtProgressBar(min=1, max=length(dt2), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  
  lin_prior <- c(brms::set_prior("normal(0, 10)", class="b"),
                 brms::set_prior("normal(0, 10)", class="Intercept"),
                 brms::set_prior("normal(0.05, 0.5)", class="sigma"))
  
  sig_prior <- c(brms::set_prior("normal(0, 10)", nlpar = "A"),
                 brms::set_prior("normal(0, 10)", nlpar = "K", lb=0.01),
                 brms::set_prior("normal(0, 10)", nlpar = "M"), 
                 brms::set_prior("normal(0, 10)", nlpar = "B"),
                 brms::set_prior("normal(0.05, 0.5)", class="sigma"))
  
  dt2 <- foreach::foreach(x = dt2,.options.snow=opts,.packages=c("data.table","brms","loo"),.export = c("lin_prior","sig_prior")) %dopar% {
    
    lin <- tryCatch(rstan::stan(file = "/Users/ul20791/Desktop/Academia/PhD/Repositories/MEWS/multivariate-glv-ews/Code/R/lin_code.stan", 
                       data = make_standata(density ~ pressure, data = x,
                                            family = gaussian(),prior = lin_prior),
                       chains = 2,iter=2000,
                       cores = 2,
                       refresh = 0,
                       control=list(adapt_delta=0.975, 
                                    max_treedepth=15)),
    error=function(e){
      return(data.frame(failed = TRUE))
    })
    
    lin2 <- brms::brm(density ~ pressure, x, empty = TRUE)
    lin2$fit <- lin
    lin2 <- rename_pars(lin2)
    lin.loo <-  ifelse(inherits(lin2,"brmsfit"),
                       loo::loo(lin2)$estimates[3,1],Inf)
    
    rm(lin,lin2)
   
    sig <- tryCatch(rstan::stan(file = "/Users/ul20791/Desktop/Academia/PhD/Repositories/MEWS/multivariate-glv-ews/Code/R/sig_code.stan", 
                                data = make_standata(brms::bf(density~A+((K-A)/(1+ exp(-1*B*(pressure-M)))),
                                                              # Nonlinear variables
                                                              A + K + B + M ~ 1,
                                                              # Nonlinear fit
                                                              nl = TRUE), data = x,
                                                     family = gaussian(),prior = sig_prior),
                                chains = 2,
                                iter=5000,
                                cores = 2,
                                refresh = 0,
                                #silent =2,
                                control=list(adapt_delta=0.975, 
                                             max_treedepth=15)),
    error=function(e){
      return(data.frame(failed = TRUE))
    })
    
    sig2 <- brms::brm(brms::bf(density~A+((K-A)/(1+ exp(-1*B*(pressure-M)))),
                               # Nonlinear variables
                               A + K + B + M ~ 1,
                               # Nonlinear fit
                               nl = TRUE), x, empty = TRUE,prior = sig_prior)
    sig2$fit <- sig
    sig2 <- rename_pars(sig2)
    sig.loo <-  ifelse(inherits(sig2,"brmsfit"),
                       loo::loo(sig2)$estimates[3,1],Inf)
    max_rhat_sig <- suppressWarnings(max(as.numeric(summary(sig2)$fixed$Rhat)))
    rm(sig,sig2)
    
    out <- data.frame(var=x$sim_id[1],
                      "model" = ifelse(lin.loo <= sig.loo & max_rhat_sig >1.05,"linear","sigmoidal")) 
    names(out) <- c(var,"model")
    
    # if(i %% 1000){
    #   Sys.sleep(30)
    #   }
    
    return(out)
    
  }
  #parallel::stopCluster(cl)
  
  close(pb)
  snow::stopCluster(cl)
  
  return(data.table::rbindlist(dt2))
}


brms_shape3 <- function(dt,var,n_cores = 4){
  
  require(foreach)
  require(brms)
  require(doFuture)
  require(progressr)
  require(inflection)
  
  doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  #cl <- parallel::makeCluster(floor(n_cores/2), type = "FORK")
  future::plan(multicore, workers = floor(n_cores/2))     ## forked parallel processing (via 'parallel')

  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  #progressr::handlers(handler_txtprogressbar(enable = TRUE))
  dt2 <- split(dt,by = var)
  
  #cl <- parallel::makePSOCKcluster(floor(n_cores/2),outfile="Log.txt",.errorhandling = "pass") # number of cores to use
  
  #cl <- parallel::makeCluster(floor(n_cores/2), type = "FORK")

 # doSNOW::registerDoSNOW(cl)
  
  pb <- progressr::progressor(steps = length(dt2))
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress=progress)
  
  
  lin_prior <- c(brms::set_prior("normal(0, 10)", class="b"),
                 brms::set_prior("normal(0, 10)", class="Intercept"),
                 brms::set_prior("normal(0.05, 0.5)", class="sigma"))
  
  sig_prior <- c(brms::set_prior("normal(0, 10)", nlpar = "A"),
                 brms::set_prior("normal(0, 10)", nlpar = "K", lb=0.01),
                 brms::set_prior("normal(0, 10)", nlpar = "M"), 
                 brms::set_prior("normal(0, 10)", nlpar = "B"),
                 brms::set_prior("normal(0.05, 0.5)", class="sigma"))
  
  
  lin <- tryCatch(brms::brm(density ~ pressure, data = dt2[[1]],
                            iter = 2000,family = gaussian(), 
                            prior = lin_prior,
                            silent =2, chains = 2,
                            thin=2,refresh=0,
                            seed = 12345, 
                            cores = 2,
                            backend = "rstan",
                            #open_progress = FALSE,
                            control=list(adapt_delta=0.975, 
                                         max_treedepth=15)),error=function(e){
                                           return(data.frame(failed = TRUE))
                                         })
  
  
  sig <-  suppressMessages(brms::brm(brms::bf(density~A+((K-A)/(1+ exp(-1*B*(pressure-M)))),
                                                       # Nonlinear variables
                                                       A + K + B + M ~ 1,
                                                       # Nonlinear fit
                                                       nl = TRUE),
                                              data = dt2[[1]],
                                              iter = 5000,
                                              prior = sig_prior,
                                              silent =2,
                                              refresh=0,
                                              chains =2,
                                              family = gaussian(),
                                              #open_progress = FALSE,
                                              thin = 2,
                                              backend = "rstan",
                                              seed = 12345, 
                                              cores = 2,
                                              control=list(adapt_delta=0.975, 
                                                           max_treedepth=15)))
  
dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","brms","loo","inflection"),.export = c("lin","sig")) %dopar% {
    
    lin.tmp <- suppressMessages(tryCatch(update(lin,newdata = x,recompile = FALSE,seed = 12345),error=function(e){
      return(data.frame(failed = TRUE))
    }))
    lin.loo <-  ifelse(inherits(lin.tmp,"brmsfit"),
                       loo::loo(lin.tmp)$estimates[3,1],Inf)
    
    rm(lin.tmp)
    
    sig.tmp <-  suppressMessages(tryCatch(update(sig,newdata = x,recompile = FALSE,seed = 12345),error=function(e){
      return(data.frame(failed = TRUE))
    }))
    sig.loo <-  ifelse(inherits(sig.tmp,"brmsfit"),
                       tryCatch(loo::loo(sig.tmp)$estimates[3,1],error = function(e){return(Inf)}),Inf)
    max_rhat_sig <- tryCatch(suppressWarnings(max(as.numeric(summary(sig.tmp)$fixed$Rhat))),
                             error=function(err){return(Inf)})
    rm(sig.tmp)
    gc()
    out <- data.frame(var=x$sim_id[1],
                      "model" = ifelse(lin.loo <= sig.loo | max_rhat_sig >1.05,"linear","sigmoidal")) 
    names(out) <- c(var,"model")
    
    out <- data.frame(var=x$sim_id[1],
                      "model" = ifelse(lin.loo <= sig.loo | max_rhat_sig >1.05,"linear","sigmoidal"),
                      "inflection_pt_incr" = inflection::bese(x[["time"]],x[["density"]],0)$iplast,
                      "inflection_pt_decr" = inflection::bese(x[["time"]],x[["density"]],1)$iplast) 
    names(out) <- c(var,"model","inflection_pt_incr","inflection_pt_decr")
    
    # if(i %% 1000){
    #   Sys.sleep(30)
    #   }
    pb()
    return(out)
    
  }
  #parallel::stopCluster(cl)
  
  #close(pb)
  #snow::stopCluster(cl)
  
  return(data.table::rbindlist(dt2))
}

INLA_shape2 <- function(dt,var,n_cores = 4){
  
  require(foreach)
  require(doFuture)
  require(progressr)
  require(inflection)
  require(doRNG)
  
  doFuture::registerDoFuture()  ## %dopar% parallelizes via future
  future::plan(multisession, workers = n_cores)     ## forked parallel processing (via 'parallel')
  
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  
  dt2 <- split(dt,by = var)
  pb <- progressr::progressor(steps = length(dt2))
  
  dt2 <- foreach::foreach(x = dt2,.packages=c("data.table","INLA","inflection")) %dorng% {
    
    lin <- INLA::inla(density ~ pressure, family="gaussian", 
                      data = x,
                      silent=2,
                      control.family = list(
                        hyper = list(
                          prec = list(
                            initial = 1,
                            fixed = F))),
                      control.compute = list(waic = TRUE))
    
    sig <- suppressMessages(tryCatch(INLA::inla(density ~ f(pressure, model="sigm"),
                                                data =  x,
                                                family = "gaussian",
                                                safe=TRUE,
                                                silent=2,
                                                control.compute = list(waic = TRUE),
                                                control.family = list(
                                                  hyper = list(
                                                    prec = list(
                                                      #initial = log(1/0.001^2),
                                                      initial = 1,
                                                      fixed = F)))),error=function(e){return(list("waic" =list("waic"=Inf)))}))
    
    inv_sig <- suppressMessages(tryCatch(INLA::inla(density ~ f(pressure, model="revsigm"),
                                                    data = x,
                                                    family = "gaussian",
                                                    safe=TRUE,
                                                    silent=2,
                                                    control.compute = list(waic = TRUE),
                                                    control.family = list(
                                                      hyper = list(
                                                        prec = list(
                                                          #initial = log(1/0.001^2),
                                                          initial = 1,
                                                          fixed = F)))),error=function(e){return(list("waic" =list("waic"=Inf)))}))
    
    out <- data.frame(var=x$sim_id[1],
                      "model" = ifelse(lin$waic$waic < sig$waic$waic | lin$waic$waic < inv_sig$waic$waic,"linear","sigmoidal"),
                      "inflection_pt_incr" = inflection::bese(x[["time"]],x[["density"]],0)$iplast,
                      "inflection_pt_decr" = inflection::bese(x[["time"]],x[["density"]],1)$iplast) 
    names(out) <- c(var,"model","inflection_pt_incr","inflection_pt_decr")
    pb()
    
    return(out)
  }
  
  return(data.table::rbindlist(dt2))
}
