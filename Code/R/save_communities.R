#############################################################################
## Preamble
#############################################################################

require(pbmcapply) #parallelisation lapply (necessary as interaction matrix generation is by trial and error)

load("Data/networks/template_motifs/adjacent_matrix_list_minimal.Rdata") #
#generic interaction matrix which can be extrapolated to higher species number
adjacent_matrix_key_minimal <- read.csv("Data/networks/template_motifs/adjacent_matrix_key_minimal.csv")
#reference key for adjacent_matrix_list_minimal to reference against

source("Code/R/A_matrix_fns.R")

###################################
## Simulate 5 species communities
###################################
pbmcapply::pbmclapply(1:10,FUN = function(i){
  #i=9
  if(max(adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,]$tlvl) ==3){
    new_spp_ls <- proportion_trophic(5,2/5,2/5,1/5,0.0)
  }else{
    new_spp_ls <- proportion_trophic(5,2/5,1/5,1/5,1/5)
  }
  
  init_A <-generalise_A(adjacent_matrix_list_minimal,i,adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,-1],new_spp_ls)
  out_file <-generate_foodwebs_small(i,init_A,new_spp_ls,300,minimal =T)
  save(out_file,file = paste("Data/networks/5_spp/web_",i,".RData",sep = ""))
  
},mc.cores = 5)

###################################
## Simulate 10 species communities
###################################
pbmcapply::pbmclapply(1:10,FUN = function(i){
  #i=9
  if(max(adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,]$tlvl) ==3){
    new_spp_ls <- proportion_trophic(10,6.01/10,3/10,0.99/10,0.0) #due R sum errors, need to add slight variability
  }else{
    new_spp_ls <- proportion_trophic(10,4.01/10,3/10,2/10,0.99/10)
  }
  
  init_A <-generalise_A(adjacent_matrix_list_minimal,i,adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,-1],new_spp_ls)
  out_file <-generate_foodwebs_small(i,init_A,new_spp_ls,300,minimal =T)
  save(out_file,file = paste("Data/networks/10_spp/web_",i,".RData",sep = ""))
  
},mc.cores = 5)

###################################
## Simulate 15 species communities
###################################
pbmcapply::pbmclapply(1:10,FUN = function(i){
  #i=9
  if(max(adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,]$tlvl) ==3){
    new_spp_ls <- proportion_trophic(15,8/15,5/15,2/15,0.0)
  }else{
    new_spp_ls <- proportion_trophic(15,6/15,4/15,3/15,2/15)
  }
  
  init_A <-generalise_A(adjacent_matrix_list_minimal,i,adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,-1],new_spp_ls)
  out_file <-generate_foodwebs_large(i,init_A,new_spp_ls,300,minimal =T)
  save(out_file,file = paste("Data/networks/15_spp/web_",i,".RData",sep = ""))
  
},mc.cores = 5)

###################################
## Simulate 20 species communities
###################################
pbmcapply::pbmclapply(1:10,FUN = function(i){
  #i=9
  if(max(adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,]$tlvl) ==3){
    new_spp_ls <- proportion_trophic(20,0.5,0.35,0.15,0.0)
  }else{
    new_spp_ls <- proportion_trophic(20,0.5,0.30,0.15,0.05)
  }
  
  init_A <-generalise_A(adjacent_matrix_list_minimal,i,adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,-1],new_spp_ls)
  out_file <-generate_foodwebs_large(i,init_A,new_spp_ls,300,minimal =T)
  save(out_file,file = paste("Data/networks/20_spp/web_",i,".RData",sep = ""))
  
},mc.cores = 5)

###################################
## Simulate 25 species communities (assuming pyramid of richnesses - decreasing up trophic levels)
###################################
pbmcapply::pbmclapply(1:10,FUN = function(i){
  if(max(adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,]$tlvl) ==3){
    new_spp_ls <- proportion_trophic(25,0.5,0.35,0.15,0.0)
  }else{
    new_spp_ls <- proportion_trophic(25,0.5,0.30,0.15,0.05)
  }
  
  init_A <-generalise_A(adjacent_matrix_list_minimal,i,adjacent_matrix_key_minimal[adjacent_matrix_key_minimal$motif == i,-1],new_spp_ls)
  out_file <-generate_foodwebs_large(i,init_A,new_spp_ls,300,minimal =T)
  save(out_file,file = paste("Data/networks/25_spp/web_",i,".RData",sep = ""))
  
},mc.cores = 5)
