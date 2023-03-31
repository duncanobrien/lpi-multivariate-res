########################################################################
## Preamble ##
########################################################################
require(EWSmethods)
require(data.table)
require(magrittr)
require(dplyr)
require(inflection)
source("Code/R/parallel_resilience_fns.R")
source("Code/R/julia_wrangle_function.R")

################################################
# Generate resilience data
################################################

motif = 1
ts_len = 70
num_spp = 5
model = "invasive"
search_effort = 0.95

resilience_models <- lapply(c("harvest","invasive"),FUN = function(model){

  stress_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/stress/",sep="")
  load_stress_files <- fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "motif_.*\\.csv$"))

  unstressed_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/unstressed/",sep="")
  load_unstressed_files <- fs::dir_ls(path = unstressed_path, regexp = paste0(unstressed_path, "motif_.*\\.csv$"))

  raw_stress_data <- as.data.table(vroom::vroom(load_stress_files[motif],show_col_types = FALSE)) %>%
    dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column"))

  hypo_trans_data <- read.csv(fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "jacobian.*\\.csv$"))) %>%
    as.data.table() %>%
    .[,community_id := paste(motif,community,sep = "_")] %>%
    setnames(.,"collapse","jac_collapse")

  raw_unstressed_data <- as.data.table(vroom::vroom(load_unstressed_files[motif],show_col_types = FALSE)) %>%
    dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column"))

  inflection_detection_stress <- julia_wrangle_dt(data = copy(raw_stress_data),
                                           n_spp = num_spp,
                                           sample = NA,
                                           motif = T) %>%
    as.data.table() %>%
    .[time %in% c((max(time)-100):max(time)) &  species == "pca1",] %>%
    .[,pressure := seq(0,unique(max_stress),by = unique(max_stress)/100),by="sim_id"] %>%
    .[,inflection_pt_incr := inflection::bese(time,density,0)$iplast,by="sim_id"] %>%
    .[,inflection_pt_decr := inflection::bese(time,density,1)$iplast,by="sim_id"] %>%
    .[,c("sim_id","inflection_pt_incr","inflection_pt_decr"),with = FALSE] %>%
    unique(.)

  resilience_stress_csv <- copy(raw_stress_data) %>%
    setnames(.,"target_spp","target_node") %>%
    .[,community_id := paste(motif,community,sep = "_")] %>%
    .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
    .[,time := seq_along(spp_1),by = "sim_id"] %>%
    merge(hypo_trans_data[,c("community_id","motif","community","jac_collapse")], 
          by =c("community_id","motif","community")) %>%
    merge(inflection_detection_stress,by="sim_id")
  
  resilience_unstressed_csv <- copy(raw_stress_data) %>%
    setnames(.,"target_spp","target_node") %>%
    .[,target_node := NA] %>% #as this node experiences no stress, set to NA
    .[,jac_collapse := NA] %>% #as above
    .[,community_id := paste(motif,community,sep = "_")] %>%
    .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
    .[,time := seq_along(spp_1),by = "sim_id"] %>%
    .[,inflection_pt_incr := max(time)] %>% #if no stress, set to last time point for downstreaming processing
    .[,inflection_pt_decr :=  max(time)] 
   
  # dt1 <- subset(resilience_stress_csv,sim_id == "1_2_1") %>%
  #   .[,grep("^spp_",names(.)) := lapply(.SD, function(x){round(x*10)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
  #   .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
  #   .[time %in% (round(unique(inflection_pt_incr))-ts_len+1):round(unique(inflection_pt_incr)),]
  #   #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
  #   
  # out_dt1 <-  merge(copy(dt1),parallel_multiJI(copy(dt1), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
  #                   by =c("time","sim_id")) %>%
  #   merge(.,parallel_uniJI(copy(dt1), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
  #         by =c("time","sim_id")) %>%
  #   merge(.,parallel_FI(copy(dt1), var = "sim_id",n_cores = 4,winsize = 25),
  #         by =c("time","sim_id")) %>%
  #   merge(.,parallel_mvi(copy(dt1), var = "sim_id",n_cores = 4,winsize = 25),
  #         by =c("time","sim_id")) 
  # 
  # dt2 <- subset(resilience_stress_csv,sim_id == "1_2_1") %>%
  #   .[,grep("^spp_",names(.)) := lapply(.SD, function(x){round(x*10)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
  #   .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
  #   .[time %in% (max(unique(time))-ts_len+1):round(unique(inflection_pt_decr)),] #filter to last time point of unstressed data and period of ts_len prior
  #   #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
  #  
  # out_dt2 <-  merge(copy(dt2),parallel_multiJI(copy(dt2), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
  #                   by =c("time","sim_id")) %>%
  #   merge(.,parallel_uniJI(copy(dt2), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
  #         by =c("time","sim_id")) %>%
  #   merge(.,parallel_FI(copy(dt2), var = "sim_id",n_cores = 4,winsize = 25),
  #         by =c("time","sim_id")) %>%
  #   merge(.,parallel_mvi(copy(dt2), var = "sim_id",n_cores = 4,winsize = 25),
  #         by =c("time","sim_id")) 
  
  if(model == "harvest"){
    resilience_stress_csv <-  resilience_stress_csv %>%
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){round(x*100)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
      .[time %in% (round(unique(inflection_pt_decr))-ts_len+1):round(unique(inflection_pt_decr)),] #filter to transition time point and period of ts_len prior
      #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
      
      out_stress_csv <- merge(copy(resilience_stress_csv),parallel_multiJI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
                              by =c("time","sim_id")) %>%
        merge(.,parallel_uniJI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
              by =c("time","sim_id")) %>%
        merge(.,parallel_FI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25),
              by =c("time","sim_id")) %>%
        merge(.,parallel_mvi(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25),
              by =c("time","sim_id")) 
    
    resilience_unstressed_csv <-  resilience_unstressed_csv %>%
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){round(x*100)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
      .[time %in% (max(unique(time))-ts_len+1):round(unique(inflection_pt_decr)),] #filter to last time point of unstressed data and period of ts_len prior
      #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
      
    out_unstressed_csv <- merge(copy(resilience_unstressed_csv),parallel_multiJI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
                            by =c("time","sim_id")) %>%
      merge(.,parallel_uniJI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
            by =c("time","sim_id")) %>%
      merge(.,parallel_FI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25),
            by =c("time","sim_id")) %>%
      merge(.,parallel_mvi(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25),
            by =c("time","sim_id")) 
    
  }else if(model == "invasive"){
    resilience_stress_csv <-  resilience_stress_csv %>%
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){round(x*10)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
      .[time %in% (round(unique(inflection_pt_incr))-ts_len+1):round(unique(inflection_pt_incr)),] 
      #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
      
    out_stress_csv <- merge(copy(resilience_stress_csv),parallel_multiJI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
                            by =c("time","sim_id")) %>%
      merge(.,parallel_uniJI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
            by =c("time","sim_id")) %>%
      merge(.,parallel_FI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25),
            by =c("time","sim_id")) %>%
      merge(.,parallel_mvi(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25),
            by =c("time","sim_id")) 
    
    resilience_unstressed_csv <-  resilience_unstressed_csv %>%
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){round(x*10)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
      .[time %in% (max(unique(time))-ts_len+1):round(unique(inflection_pt_decr)),] #filter to last time point of unstressed data and period of ts_len prior
      #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
     
    out_unstressed_csv <- merge(copy(resilience_unstressed_csv),parallel_multiJI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
                                by =c("time","sim_id")) %>%
      merge(.,parallel_uniJI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
            by =c("time","sim_id")) %>%
      merge(.,parallel_FI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25),
            by =c("time","sim_id")) %>%
      merge(.,parallel_mvi(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25),
            by =c("time","sim_id")) 
  }
    out_stress_csv <- out_stress_csv %>%
      .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
           "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
           "multiJI","mean_uniJI","max_uniJI","FI","mvi"),with = FALSE] %>% #select columns of interest
      .[,model := model] %>% #add additional metadata
      .[,stressed := 1] %>%
      .[,ts_length := ts_len] %>%
      .[,search_effort := search_effort]
    
    out_unstressed_csv <- out_unstressed_csv %>%
      .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
           "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
           "multiJI","mean_uniJI","max_uniJI","FI","mvi"),with = FALSE] %>% #select columns of interest
      .[,model := model] %>% #add additional metadata
      .[,stressed := 0] %>%
      .[,ts_length := ts_len] %>%
      .[,search_effort := search_effort]
    
    out <- rbind(out_stress_csv,out_unstressed_csv)
    return(out)

    }) %>%
  data.table::rbindlist()