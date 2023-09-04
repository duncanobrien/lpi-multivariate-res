########################################################################
## Preamble ##
########################################################################
require(EWSmethods)
require(data.table)
require(magrittr)
require(dplyr)
require(inflection)
source("Code/R/auxiliary_functions/parallel_resilience_fns.R")
source("Code/R/auxiliary_functions/julia_wrangle_function.R")

################################################
# Generate resilience data
################################################

#motifs = c(1,2,6,9,10)
motifs = 1
#ts_len = 70
ts_len = seq(10,70,15)
#ts_len[1] = 11
#num_spp = 5
#num_spps = seq(5,25,5)
num_spps = c(5,15,25)
model = "invasive"
#search_effort = 0.95
search_effort = seq(0.1,1.0,0.1)
model = c("harvest","invasive")

parameter_space <- expand.grid("motif" = motifs,
                               "num_spp" = num_spps,
                               "ts_len" = ts_len,
                               "search_effort" = search_effort,
                               "model" = c("harvest","invasive")) %>%
  dplyr::arrange(motif,num_spp,ts_len,search_effort) %>%
  #base::split(interaction(.$motif,.$num_spp))
#base::split(interaction(.$motif,.$num_spp,.$ts_len))
base::split(interaction(.$motif,.$num_spp,.$model))

sample_id <- expand.grid("motif" = 1,
                         "community" = 1:30,
                         "sim" = 1:25) |>
  dplyr::mutate(sample_id = paste(motif,community,sim,sep = "_")) |>
  dplyr::select(sample_id) |> unlist() |> unname()


i=1
motif_data <- parameter_space[[6]]

lapply(parameter_space, function(motif_data){
  
  motif = unique(motif_data$motif)
  num_spp =  unique(motif_data$num_spp)
  #ts_len = unique(motif_data$ts_len)
  mod = as.character(unique(motif_data$model))
  
  resilience_models <- lapply(c(mod),FUN = function(model){
    
    print(paste0("Preparing ",model," data"))
    
    stress_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/stress/",sep="") #define stressed model path
    load_stress_files <- fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "motif_.*\\.csv$"))
    
    unstressed_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/unstressed/",sep="") #define unstressed model path
    load_unstressed_files <- fs::dir_ls(path = unstressed_path, regexp = paste0(unstressed_path, "motif_.*\\.csv$"))
    
    raw_stress_data <- as.data.table(vroom::vroom(load_stress_files[grepl(paste0("_",motif,"\\."),load_stress_files)],show_col_types = FALSE)) %>%
      dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column")) #load stressed models and rename species columns
    
    hypo_trans_data <- read.csv(fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "jacobian.*\\.csv$"))) %>% #load stressed models' jacobian info
      as.data.table() %>%
      .[,community_id := paste(motif,community,sep = "_")] %>% #set simulation identifier
      setnames(.,"collapse","jac_collapse") #rename for readability
    
    raw_unstressed_data <- as.data.table(vroom::vroom(load_unstressed_files[grepl(paste0("_",motif,"\\."),load_unstressed_files)],show_col_types = FALSE)) %>%
      dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column")) #load unstressed models and rename species columns
    
    inflection_detection_stress <- julia_wrangle_dt(data = copy(raw_stress_data),
                                                    n_spp = num_spp,
                                                    sample = NA,
                                                    motif = T) %>% #wrangle stressed models to generate PCA and provide simulation identifer
      as.data.table() %>%
      .[time %in% c((max(time)-100):max(time)) &  species == "pca1",] %>% #add time vector
      .[,pressure := seq(0,unique(max_stress),by = unique(max_stress)/100),by="sim_id"] %>% #and link this time vector to control parameter
      .[,inflection_pt_incr := inflection::bese(time,density,0,doparallel = F)$iplast,by="sim_id"] %>% #identify positive turnpoint 
      .[,inflection_pt_decr := inflection::bese(time,density,1,doparallel = F)$iplast,by="sim_id"] %>% #and negative turnpoint
      .[,grep("^inflection",names(.)) := lapply(.SD, function(x){ifelse(is.na(x),return(max(time)),return(x))}), .SDcols = grep("^inflection",names(.)), by = "sim_id"] %>% #if no inflection point found, set to last time point
      .[,c("sim_id","inflection_pt_incr","inflection_pt_decr"),with = FALSE] %>%
      unique(.) #filter to remove duplicates
    
    resilience_stress_csv <- copy(raw_stress_data) %>% #merge the raw data with inflection point data
      setnames(.,"target_spp","target_node") %>%
      .[,community_id := paste(motif,community,sep = "_")] %>%
      .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
      .[,time := seq_along(spp_1),by = "sim_id"] %>%
      merge(hypo_trans_data[,c("community_id","motif","community","jac_collapse")], 
            by =c("community_id","motif","community")) %>%
      merge(inflection_detection_stress,by="sim_id") %>%
      .[,first_neg := min(ifelse(apply(.SD, 1,min) < 0 , time, max(time))),.SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #due to erroneous negative values from simulation solver, must crop prior to first negative value
      .[.[,.I[time < first_neg], by = "sim_id"]$V1] %>% #filter prior to first negative value 
      .[,first_neg:=NULL] %>%
      .[sim_id %in% sample_id,]
    
    rm(inflection_detection_stress) #remove object to release memory 
    
    resilience_unstressed_csv <- copy(raw_unstressed_data) %>% #repeat the inflection point analysis but for the unstressed models
      setnames(.,"target_spp","target_node") %>%
      .[,target_node := NA] %>% #as no node experiences stress, set to NA
      .[,jac_collapse := NA] %>% #as above but for the jacobian
      .[,community_id := paste(motif,community,sep = "_")] %>%
      .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
      .[,time := seq_along(spp_1),by = "sim_id"] %>%
      .[,first_neg := min(ifelse(apply(.SD, 1,min) < 0 , time, max(time))),.SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #due to erroneous negative values from simulation solver, must crop prior to first negative value
      .[.[,.I[time < first_neg], by = "sim_id"]$V1] %>% #filter prior to first negative value 
      .[,first_neg:=NULL] %>% 
      .[,inflection_pt_incr := max(time), by = "sim_id"] %>% #as no stress, set inflection point to last time point for downstreaming processing
      .[,inflection_pt_decr :=  max(time), by = "sim_id"] %>%
      .[sim_id %in% sample_id,]
    
    rm(raw_stress_data,raw_unstressed_data) #remove object to release memory 
    
    inflec_col <- ifelse(model =="harvest","inflection_pt_decr","inflection_pt_incr") #harvesting model expected to negatively inflect, the invasive model, positively
    
    save_resilience_data <- lapply(seq_len(nrow(motif_data)), function(df){
      #save_resilience_data <- lapply(1:2, function(df){
      closeAllConnections()
      parameter_data <- motif_data[df,]
      ts_len = parameter_data$ts_len
      search_effort = parameter_data$search_effort
      print(paste(motif,num_spp,model,ts_len,search_effort,sep="_"))
      
      winsize = ifelse(ts_len == 10, 60,50)
      
      tmp_stress_csv <-  copy(resilience_stress_csv) %>%
        .[,grep("^spp_",names(.)) := lapply(.SD, function(x){ifelse(model =="harvest",return(as.integer(x*1000)),return(as.integer(x*10)))}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
        .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
        .[.[, .I[time %in% (round(unique(get(inflec_col)))-ts_len+1):round(unique(get(inflec_col)))], by = "sim_id"]$V1] #conditionally filter (by ts_len) by simulation identifier
      
      out_stress_csv <- merge(copy(tmp_stress_csv),parallel_mvi(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize), #calculate each resilience metric for the stressed models
                              by =c("time","sim_id")) %>% 
        merge(.,parallel_multiJI(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,sample_spp = TRUE, winsize = winsize),
                              by =c("time","sim_id"),all.x = TRUE) %>% 
        merge(.,parallel_uniJI(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize,E = 1,tau = round(ts_len*-0.1)),
              by =c("time","sim_id"),all.x = TRUE) %>%
        merge(.,parallel_FI(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize,TL=75),
              by =c("time","sim_id"),all.x = TRUE) %>%
        merge(.,parallel_multiAR(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize),
              by =c("time","sim_id"),all.x = TRUE) %>%
        .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
             "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
             "multiJI",names(.)[grep("^uniJI_",names(.))],"mean_uniJI","max_uniJI","FI","mvi","multiAR"),with = FALSE] %>% #select columns of interest
        .[,model := model] %>% #add additional metadata
        .[,stressed := 1] %>%
        .[,ts_length := ts_len] %>%
        .[,search_effort := search_effort] %>%
        .[order(sim_id,time),] 
      
      tmp_unstressed_csv <-  copy(resilience_unstressed_csv) %>%
        .[,grep("^spp_",names(.)) := lapply(.SD, function(x){ifelse(model =="harvest",return(as.integer(x*1000)),return(as.integer(x*10)))}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
        .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
        .[.[, .I[time %in% (round(unique(get(inflec_col)))-ts_len+1):round(unique(get(inflec_col)))], by = "sim_id"]$V1] #conditionally filter (by ts_len) by simulation identifier
      #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
      
      out_unstressed_csv <-  merge(copy(tmp_unstressed_csv),parallel_mvi(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize), #calculate each resilience metric for the unstressed models
                                   by =c("time","sim_id")) %>%
        merge(.,parallel_multiJI(copy(tmp_unstressed_csv), var = "sim_id",n_cores =9,sample_spp = TRUE,winsize = winsize),
                                  by =c("time","sim_id"),all.x = TRUE) %>% 
        merge(.,parallel_uniJI(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize,E = 1,tau = round(ts_len*-0.1)),
              by =c("time","sim_id"),all.x = TRUE) %>%
        merge(.,parallel_FI(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize),
              by =c("time","sim_id"),all.x = TRUE) %>%
        merge(.,parallel_multiAR(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize),
              by =c("time","sim_id"),all.x = TRUE) %>%
        .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
             "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
             "multiJI",names(.)[grep("^uniJI_",names(.))],"mean_uniJI","max_uniJI","FI","mvi","multiAR"),with = FALSE] %>% #select columns of interest
        .[,model := model] %>% #add additional metadata
        .[,stressed := 0] %>%
        .[,ts_length := ts_len] %>%
        .[,search_effort := search_effort] %>%
        .[order(sim_id,time),] 
      
      out <- rbind(copy(out_stress_csv),copy(out_unstressed_csv)) #merge stressed and unstressed data for saving
      rm(out_stress_csv,out_unstressed_csv,tmp_stress_csv,tmp_unstressed_csv)
      gc()
      return(out)
    }) %>%
      data.table::rbindlist() 
    
    return(save_resilience_data)
  }) %>%
    data.table::rbindlist() %>%
    .[,motif := motif] %>%
    .[,n_spp := num_spp]
  
  resilience_summary_data <- resilience_models %>%
    melt(.,measure.vars = c("multiJI",names(.)[grep("^uniJI_",names(.))],"mean_uniJI","max_uniJI","FI","mvi","multiAR"),
         variable.name = "metric",value.name = "metric_value") %>% #pivot resilience metrics longer
    .[,trend :=  tryCatch(cor.test(time, metric_value, alternative = c("two.sided"), method = c("kendall"), conf.level = 0.95,na.action = na.omit)$estimate,error = function(err){as.numeric(NA)}),
      by = c("ts_length","search_effort","sim_id","metric", "stressed","model")] %>% #estimate Kendall tau correlation through time
    .[,corrected_trend :=  ifelse(metric == "FI",trend*-1,trend), #as FI is expected to decrease, use it's absolute value to be comparable
      by = c("ts_length","search_effort","sim_id","metric", "stressed","model")] %>%
    .[,threshold_crossed :=  ifelse(metric %in% c("FI","mvi") | all(is.na(metric_value)), as.numeric(NA),
                                    ifelse(metric %in% "multiAR" & any(metric_value > 0), 1,
                                    ifelse(any(metric_value > 1), 1,0))),by = c("ts_length","search_effort","sim_id","metric", "stressed","model")] %>% #for the Jacobian indices, a value >1 = unstable
    .[,.SD[1],  c("ts_length","search_effort","sim_id","metric", "stressed","model")] %>% #select unique
    .[,c("sim_id","metric","trend","corrected_trend","threshold_crossed","target_node","jac_collapse","model","stressed","ts_length","search_effort","motif","n_spp")] %>% #drop time
    .[order(model,stressed,sim_id),] #order for readability
  
  #write.csv(save_resilience_data,file = paste0("Data/resilience/full/motif",motif,"_",num_spp,"spp.csv"))
  #write.csv(resilience_summary_data,file = paste0("Data/resilience/summary/motif",motif,"_",num_spp,"spp.csv"))
  data.table::fwrite(resilience_models,file = paste0("Data/resilience/full/motif",motif,"_",num_spp,"_",mod,".csv"),compress = "gzip",nThread = 8)
  saveRDS(resilience_models,file = paste0("Data/resilience/full/motif",motif,"_",num_spp,"_",mod,".rds"))
  
  data.table::fwrite(resilience_summary_data,file = paste0("Data/resilience/summary/summary_motif",motif,"_",num_spp,"_",mod,".csv"),compress = "gzip",nThread = 8)
  saveRDS(resilience_summary_data,file = paste0("Data/resilience/summary/summary_motif",motif,"_",num_spp,"_",mod,".rds"))
  
  rm(resilience_summary_data,resilience_models)
  gc()
  
})
