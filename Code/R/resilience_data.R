########################################################################
## Preamble ##
########################################################################
#require(EWSmethods)
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

  stress_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/stress/",sep="") #define stressed model path
  load_stress_files <- fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "motif_.*\\.csv$"))

  unstressed_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/unstressed/",sep="") #define unstressed model path
  load_unstressed_files <- fs::dir_ls(path = unstressed_path, regexp = paste0(unstressed_path, "motif_.*\\.csv$"))

  raw_stress_data <- as.data.table(vroom::vroom(load_stress_files[motif],show_col_types = FALSE)) %>%
    dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column")) #load stressed models and rename species columns

  hypo_trans_data <- read.csv(fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "jacobian.*\\.csv$"))) %>% #load stressed models' jacobian info
    as.data.table() %>%
    .[,community_id := paste(motif,community,sep = "_")] %>% #set simulation identifier
    setnames(.,"collapse","jac_collapse") #rename for readability

  raw_unstressed_data <- as.data.table(vroom::vroom(load_unstressed_files[motif],show_col_types = FALSE)) %>%
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
    merge(inflection_detection_stress,by="sim_id")
  
  rm(inflection_detection_stress) #remove object to release memory 
  
  resilience_unstressed_csv <- copy(raw_unstressed_data) %>% #repeat the inflection point analysis but for the unstressed models
    setnames(.,"target_spp","target_node") %>%
    .[,target_node := NA] %>% #as no node experiences stress, set to NA
    .[,jac_collapse := NA] %>% #as above but for the jacobian
    .[,community_id := paste(motif,community,sep = "_")] %>%
    .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
    .[,time := seq_along(spp_1),by = "sim_id"] %>%
    .[,inflection_pt_incr := max(time)] %>% #as no stress, set inflection point to last time point for downstreaming processing
    .[,inflection_pt_decr :=  max(time)] 
   
  rm(raw_stress_data,raw_unstressed_data) #remove object to release memory 
  
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
  
  inflec_col <- ifelse(model =="harvest","inflection_pt_decr","inflection_pt_incr") #harvesting model expected to positively inflect, the invasive model, negatively
  
  resilience_stress_csv <-  resilience_stress_csv %>%
    .[,grep("^spp_",names(.)) := lapply(.SD, function(x){ifelse(model =="harvest",return(round(x*100)),return(round(x*10)))}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
    .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
    .[.[, .I[time %in% (round(unique(get(inflec_col)))-ts_len+1):round(unique(get(inflec_col)))], by = "sim_id"]$V1] #conditionally filter (by ts_len) by simulation identifier
  #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
  
  # sample_id <- unique(resilience_stress_csv$sim_id)
  # tt <- copy(resilience_stress_csv) %>%
  #   subset(.,sim_id %in% sample_id[1:100])
  # out_stress_csv <- merge(copy(tt),parallel_multiJI(copy(tt), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
  #                         by =c("sim_id","time")) %>%
  #     merge(.,parallel_uniJI(copy(tt), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
  #           by =c("time","sim_id")) %>%
  #     merge(.,parallel_FI(copy(tt), var = "sim_id",n_cores = 4,winsize = 25),
  #           by =c("time","sim_id")) %>%
  #     merge(.,parallel_mvi(copy(tt), var = "sim_id",n_cores = 4,winsize = 25),
  #           by =c("time","sim_id"))
  # 
  # kk <- copy(resilience_unstressed_csv) %>%
  #   subset(.,sim_id %in% sample_id[1:100])
  # out_unstressed_csv <- merge(copy(kk),parallel_multiJI(copy(kk), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
  #                         by =c("sim_id","time")) %>%
  #   merge(.,parallel_uniJI(copy(kk), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
  #         by =c("time","sim_id")) %>%
  #   merge(.,parallel_FI(copy(kk), var = "sim_id",n_cores = 4,winsize = 25),
  #         by =c("time","sim_id")) %>%
  #   merge(.,parallel_mvi(copy(kk), var = "sim_id",n_cores = 4,winsize = 25),
  #         by =c("time","sim_id"))

  # out_stress_csv <- merge(copy(tt),parallel_multiJI(copy(tt), var = "sim_id",n_cores = 4,winsize = 25,scale = T),
  #                         by =c("sim_id","time"))
  # 
  # save(tt,file = "/Users/ul20791/Downloads/tt_debug.RData")
  # load(file = "/Users/ul20791/Downloads/tt_debug.RData")
  
  out_stress_csv <- merge(copy(resilience_stress_csv),parallel_multiJI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
                          by =c("time","sim_id")) %>% #calculate each resilience metric for the stressed models
    merge(.,parallel_uniJI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
          by =c("time","sim_id")) %>%
    merge(.,parallel_FI(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25),
          by =c("time","sim_id")) %>%
    merge(.,parallel_mvi(copy(resilience_stress_csv), var = "sim_id",n_cores = 4,winsize = 25),
          by =c("time","sim_id")) 
  
resilience_unstressed_csv <-  resilience_unstressed_csv %>%
    .[,grep("^spp_",names(.)) := lapply(.SD, function(x){ifelse(model =="harvest",return(round(x*100)),return(round(x*10)))}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
    .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
    .[.[, .I[time %in% (round(unique(get(inflec_col)))-ts_len+1):round(unique(get(inflec_col)))], by = "sim_id"]$V1] #conditionally filter (by ts_len) by simulation identifier
  #.[time >= max(time)-ts_len+1,] %>% #stress begins at t100 but we crop here to representative lengths of LPI
  
  out_unstressed_csv <- merge(copy(resilience_unstressed_csv),parallel_multiJI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE),
                              by =c("time","sim_id")) %>% #calculate each resilience metric for the stressed models
    merge(.,parallel_uniJI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25,scale = TRUE,E = 1),
          by =c("time","sim_id")) %>%
    merge(.,parallel_FI(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25),
          by =c("time","sim_id")) %>%
    merge(.,parallel_mvi(copy(resilience_unstressed_csv), var = "sim_id",n_cores = 4,winsize = 25),
          by =c("time","sim_id")) 
  
    out_stress_csv <- out_stress_csv %>%
      .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
           "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
           "multiJI","mean_uniJI","max_uniJI","FI","mvi"),with = FALSE] %>% #select columns of interest
      .[,model := model] %>% #add additional metadata
      .[,stressed := 1] %>%
      .[,ts_length := ts_len] %>%
      .[,search_effort := search_effort] %>%
      .[order(sim_id,time),] 
    
    out_unstressed_csv <- out_unstressed_csv %>%
      .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
           "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
           "multiJI","mean_uniJI","max_uniJI","FI","mvi"),with = FALSE] %>% #select columns of interest
      .[,model := model] %>% #add additional metadata
      .[,stressed := 0] %>%
      .[,ts_length := ts_len] %>%
      .[,search_effort := search_effort] %>%
      .[order(sim_id,time),] 

    
    out <- rbind(out_stress_csv,out_unstressed_csv) #merge stressed and unstressed data for saving
    return(out)

    }) %>%
  data.table::rbindlist()

resilience_summary_data <- copy(out) %>%
  melt(.,measure.vars = c("multiJI","mean_uniJI","max_uniJI","FI","mvi"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot resilience metrics longer
  .[,trend :=  cor.test(time, metric_value, alternative = c("two.sided"), method = c("kendall"), conf.level = 0.95,na.action = na.omit)$estimate,
    by = c("sim_id","metric", "model","stressed")] %>% #estimate Kendall tau correlation through time
  .[,corrected_trend :=  ifelse(metric == "FI",abs(trend),trend), #as FI is expected to decrease, use it's absolute value to be comparable
    by = c("sim_id","metric", "model","stressed")] %>%
  .[,threshold_crossed :=  ifelse(metric %in% c("FI","mvi"), NA,
      ifelse(any(metric_value > 1), 1,0)),by = c("sim_id","metric", "model","stressed")] %>% #for the Jacobian indices, a value >1 = unstable
  .[,.SD[1],  c("sim_id","metric", "model","stressed")] %>% #drop time
  .[,c("sim_id","metric","trend","corrected_trend","threshold_crossed","target_node","jac_collapse","model","stressed","ts_length","search_effort")] %>%
  .[order(model,stressed,sim_id),] #order for readability
  