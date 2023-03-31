#' Wrangle Raw Time Series Files For Shape Fitting
#'
#' Reformats simulated generalised Lotka-Volterra models outputs in to an appropriate form for INLA shape model fitting.
#'
#' @param data A dataframe of model outputs.
#' @param n_spp Numeric.  The number of species/nodes in the simulation.
#' @param sample A numeric value stating the number of simulations to sample from. If \code{NA}, all simulations per motif are kept
#' @param motif Boolean. If varying motifs are present in \code{data}, \code{motif} must equal \code{TRUE},  
#' @returns A dataframe of wrangled simulation densities and additional principal component. 
#'

julia_wrangle_dt <- function(data,n_spp,sample = 1,motif = T){
  require(data.table)
  require(magrittr)
  if(isTRUE(motif)){
    out <-  as.data.table(data) %>%
      .[,time := seq_along(spp_1),by = c("motif","community","sim")] 
    
    if(is.na(sample)){
      out <- out %>%
        .[,.SD[!any(colSums(.SD) == 0)],by=c("motif","community","sim")] %>%
        .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
        .[,pca1 :=  prcomp(.SD,scale. = T)$x[,1],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        .[,pca2 :=  prcomp(.SD,scale. = T)$x[,2],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        data.table::melt(., measure.vars = c(paste("spp",1:n_spp,sep = "_"),"pca1","pca2"),
                         variable.name = "species", value.name = "density")
    }else{
      out <- out %>%
        .[,.SD[!any(colSums(.SD) == 0)],by=c("motif","community","sim")] %>%
        .[, list(data=list(.SD)), by=c("motif","community","sim")] %>%
        .[.[, .I[sample(.N, sample)],by=c("motif","community")]$V1] %>%
        .[,sim := seq(1:sample),by=c("motif","community")]  %>%
        .[,  data.table::rbindlist(data), by = c("motif","community","sim")] %>%
        .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
        .[,pca1 :=  prcomp(.SD,scale. = T)$x[,1],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        .[,pca2 :=  prcomp(.SD,scale. = T)$x[,2],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        data.table::melt(., measure.vars = c(paste("spp",1:n_spp,sep = "_"),"pca1","pca2"),
                         variable.name = "species", value.name = "density")
      
    }
    
  }
  if(isFALSE(motif)){
    
    out <-  as.data.table(data) %>%
      .[,time := seq_along(spp_1),by = c("community","sim")]
    
    if(is.na(sample)){
      out <- out %>%
        .[,.SD[!any(colSums(.SD) == 0)],by=c("community","sim")] %>%
        .[,sim_id := paste(community,sim,sep = "_")] %>%
        .[,pca1 :=  prcomp(.SD,scale. = T)$x[,1],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        .[,pca2 :=  prcomp(.SD,scale. = T)$x[,2],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        data.table::melt(., measure.vars = c(paste("spp",1:n_spp,sep = "_"),"pca1","pca2"),
                         variable.name = "species", value.name = "density")
    }else{
      out <- out %>%
        .[,.SD[!any(colSums(.SD) == 0)],by=c("community","sim")] %>%
        .[, list(data=list(.SD)), by=c("community","sim")] %>%
        .[.[, .I[sample(.N, sample)],by=c("community")]$V1] %>%
        .[,sim := seq(1:sample),by=c("community")]  %>%
        .[,  data.table::rbindlist(data), by = c("community","sim")] %>%
        .[,sim_id := paste(community,sim,sep = "_")] %>%
        .[,pca1 :=  prcomp(.SD,scale. = T)$x[,1],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        .[,pca2 :=  prcomp(.SD,scale. = T)$x[,2],.SDcols = paste("spp",1:n_spp,sep = "_"),by = "sim_id"] %>%
        data.table::melt(., measure.vars = c(paste("spp",1:n_spp,sep = "_"),"pca1","pca2"),
                         variable.name = "species", value.name = "density")
      
    }
  }
  return(as.data.frame(out))
}
