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

ts_len = seq(10,70,15)
model = "invasive"
search_effort = seq(0.2,1.0,0.4)
model = "invasive"

node_comb <- t(combn(paste0("spp_",1:15),5)) |>
  as.data.frame() |>
  #slice_sample(n=1000) |>
  `colnames<-`(c(paste0("node",1:5)))

spp_parameter_space <- expand.grid("motif" = 1,
                               "num_spp" = 15,
                               "ts_len" = ts_len,
                               "search_effort" = search_effort,
                               "model" = "invasive",
                               stringsAsFactors = FALSE) %>%
  dplyr::arrange(motif,num_spp,ts_len,search_effort)

node_space <- merge(spp_parameter_space,node_comb) %>%
  base::split(interaction(.$ts_len,.$search_effort))

motif = 1
num_spp = 15
mod = model

sample_id <- expand.grid("motif" = 1,
                         "community" = 1:30,
                         "sim" = 1:25) |>
  dplyr::mutate(sample_id = paste(motif,community,sim,sep = "_")) |>
  dplyr::select(sample_id) |> unlist() |> unname()

spp_resilience_models <- lapply(c(mod),FUN = function(model){
  
  print(paste0("Preparing ",model," data"))
  
  stress_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/stress/",sep="") #define stressed model path
  load_stress_files <- fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "motif_.*\\.csv$"))
  
  unstressed_path <- paste("Data/simulations/",num_spp,"_spp/",model,"/unstressed/",sep="") #define unstressed model path
  load_unstressed_files <- fs::dir_ls(path = unstressed_path, regexp = paste0(unstressed_path, "motif_.*\\.csv$"))
  
  raw_stress_data <- as.data.table(vroom::vroom(load_stress_files[grepl(paste0("_",motif,"\\."),load_unstressed_files)],show_col_types = FALSE)) %>%
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
  
  inflec_col <- ifelse(model =="harvest","inflection_pt_decr","inflection_pt_incr") #harvesting model expected to positively inflect, the invasive model, negatively
  
  node_loop <- lapply(node_space, function(kk){
    all_spp_names <- paste0("spp_",1:num_spp)
    print(paste(motif,num_spp,model,ts_len,search_effort,sep="_"))
    
  save_resilience_data <- pbapply::pblapply(seq_len(nrow(kk)), FUN = function(df){
    #save_resilience_data <- lapply(1:2, function(df){
    closeAllConnections()
    parameter_data <- kk[df,]
    ts_len = parameter_data$ts_len
    search_effort = parameter_data$search_effort
    target_spp = unlist(parameter_data[,grep("^node",colnames(parameter_data))])
    non_target_spp = setdiff(all_spp_names,target_spp)  
    
    winsize = ifelse(ts_len == 10, 60,50)
    
    tmp_stress_csv <-  copy(resilience_stress_csv) %>%
      .[,c(non_target_spp) := NULL] %>%
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){ifelse(model =="harvest",return(as.integer(x*1000)),return(as.integer(x*10)))}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
      .[.[, .I[time %in% (round(unique(get(inflec_col)))-ts_len+1):round(unique(get(inflec_col)))], by = "sim_id"]$V1] #conditionally filter (by ts_len) by simulation identifier

      out_stress_csv <- merge(copy(tmp_stress_csv),parallel_multiJI(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,sample_spp = TRUE, winsize = winsize),
                            by =c("time","sim_id")) %>% #calculate each resilience metric for the stressed models
      merge(.,parallel_uniJI(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize,E = 1,tau = round(ts_len*-0.1)),
            by =c("time","sim_id")) %>%
      merge(.,parallel_FI(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize,TL=75),
            by =c("time","sim_id")) %>%
      merge(.,parallel_mvi(copy(tmp_stress_csv), var = "sim_id",n_cores = 9,winsize = winsize),
            by =c("time","sim_id")) %>%
      .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
           "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
           "multiJI",names(.)[grep("^uniJI_",names(.))],"mean_uniJI","max_uniJI","FI","mvi"),with = FALSE] %>% #select columns of interest
      .[,model := model] %>% #add additional metadata
      .[,stressed := 1] %>%
      .[,ts_length := ts_len] %>%
      .[,search_effort := search_effort] %>%
      .[,sampled_spp := paste(tstrsplit(target_spp,"spp_")[[2]],collapse = "_")] %>%
      .[order(sim_id,time),] 
    
    tmp_unstressed_csv <-  copy(resilience_unstressed_csv) %>%
      .[,c(non_target_spp) := NULL] %>%
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){ifelse(model =="harvest",return(as.integer(x*1000)),return(as.integer(x*10)))}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #convert densities to integers of appropriate scale to LPI
      .[,grep("^spp_",names(.)) := lapply(.SD, function(x){rbinom(length(x), x, search_effort)}), .SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #introduce sampling error
      .[.[, .I[time %in% (round(unique(get(inflec_col)))-ts_len+1):round(unique(get(inflec_col)))], by = "sim_id"]$V1] #conditionally filter (by ts_len) by simulation identifier

    out_unstressed_csv <- merge(copy(tmp_unstressed_csv),parallel_multiJI(copy(tmp_unstressed_csv), var = "sim_id",n_cores =9,sample_spp = TRUE,winsize = winsize),
                                by =c("time","sim_id")) %>% #calculate each resilience metric for the stressed models
      merge(.,parallel_uniJI(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize,E = 1,tau = round(ts_len*-0.1)),
            by =c("time","sim_id")) %>%
      merge(.,parallel_FI(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize),
            by =c("time","sim_id")) %>%
      merge(.,parallel_mvi(copy(tmp_unstressed_csv), var = "sim_id",n_cores = 9,winsize = winsize),
            by =c("time","sim_id")) %>%
      .[,c("sim_id",names(.)[grep("^spp_",names(.))],"target_node",
           "time","jac_collapse","inflection_pt_incr","inflection_pt_decr",
           "multiJI",names(.)[grep("^uniJI_",names(.))],"mean_uniJI","max_uniJI","FI","mvi"),with = FALSE] %>% #select columns of interest
      .[,model := model] %>% #add additional metadata
      .[,stressed := 0] %>%
      .[,ts_length := ts_len] %>%
      .[,search_effort := search_effort] %>%
      .[,sampled_spp := paste(tstrsplit(target_spp,"spp_")[[2]],collapse = "_")] %>%
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

spp_resilience_summary_data <- resilience_models %>%
  melt(.,measure.vars = c("multiJI",names(.)[grep("^uniJI_",names(.))],"mean_uniJI","max_uniJI","FI","mvi"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot resilience metrics longer
  .[,trend :=  tryCatch(cor.test(time, metric_value, alternative = c("two.sided"), method = c("kendall"), conf.level = 0.95,na.action = na.omit)$estimate,error = function(err){NA}),
    by = c("sampled_spp","ts_length","search_effort","sim_id","metric", "stressed","model")] %>% #estimate Kendall tau correlation through time
  .[,corrected_trend :=  ifelse(metric == "FI",trend*-1,trend), #as FI is expected to decrease, use it's absolute value to be comparable
    by = c("sampled_spp","ts_length","search_effort","sim_id","metric", "stressed","model")] %>%
  .[,threshold_crossed :=  ifelse(metric %in% c("FI","mvi") | all(is.na(metric_value)), NA,
                                  ifelse(any(metric_value > 1), 1,0)),by = c("sampled_spp","ts_length","search_effort","sim_id","metric", "stressed","model")] %>% #for the Jacobian indices, a value >1 = unstable
  .[,.SD[1],  c("sampled_spp","ts_length","search_effort","sim_id","metric", "stressed","model")] %>% #select unique
  .[,c("sim_id","metric","trend","corrected_trend","threshold_crossed","target_node","jac_collapse","model","stressed","sampled_spp","ts_length","search_effort","motif","n_spp")] %>% #drop time
  .[order(model,stressed,sim_id),] #order for readability

data.table::fwrite(resilience_models,file = paste0("Data/resilience/full/motif",motif,"_",num_spp,"_",mod,".csv"),compress = "gzip",nThread = 8)
saveRDS(resilience_models,file = paste0("Data/resilience/full/motif",motif,"_",num_spp,"_",mod,".rds"))

data.table::fwrite(resilience_summary_data,file = paste0("Data/resilience/summary/summary_motif",motif,"_",num_spp,"_",mod,".csv"),compress = "gzip",nThread = 8)
saveRDS(resilience_summary_data,file = paste0("Data/resilience/summary/summary_motif",motif,"_",num_spp,"_",mod,".rds"))

rm(resilience_summary_data,resilience_models)
gc()

})

