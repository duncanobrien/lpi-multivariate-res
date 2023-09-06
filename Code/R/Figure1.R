# Figure 1

require(ggplots)
require(patchwork)
require(magick)
require(igraph)
require(ggraph)
require(tidygraph)
require(data.table)
require(magrittr)
source("Code/R/auxiliary_functions/julia_wrangle_function.R")

motif = 1 
model = "invasive"
n_spp = 5
simulation_id <- "1_2_1" #example simulation

#######################
# Load in interaction matrices
#######################

matrix_path <- paste("Data/networks/",n_spp,"_spp/","web_",motif,".RData",sep="") #define stressed model path
load(matrix_path)

simulation_info <- out_file[[as.numeric(strsplit(simulation_id,split = "_")[[1]][2])]] #extract just the example simulation

interaction_matrix <- simulation_info$A_matrix

tlvl_matrix <- as.data.frame(simulation_info$tlvl)

for(j in unique(tlvl_matrix$species)){
  if(sum(tlvl_matrix$tlvl == subset(tlvl_matrix,species == j)$tlvl) > 1){
    interaction_matrix[j,which(tlvl_matrix$tlvl == subset(tlvl_matrix,species == j)$tlvl)] <- 1
  }
} #convert intra-trophic interactions to positive to ensure they are present in the plot 

interaction_matrix <- abs(ceiling(interaction_matrix)) #extract just the positive interactions
diag(interaction_matrix) = 1 #also include intraspecific interaction

graph5 <- igraph::graph_from_adjacency_matrix(interaction_matrix,
                                              mode="directed", diag=T,weighted = T) #create graph

igraph::V(graph5)$tlvl <- simulation_info$tlvl[,1] #label nodes with trophic levels
igraph::V(graph5)$species <- simulation_info$tlvl[,2] #label nodes with species levels
igraph::E(graph5)$shared <- c("self","intra","intra","self","inter","inter",
                              "inter","inter","inter","inter","inter") #label edges with interaction types

example_graph <- tidygraph::as_tbl_graph(graph5, directed = FALSE) #make tidy

#######################
# Load in raw simulations
#######################
stress_path <- paste("Data/simulations/",n_spp,"_spp/",model,"/stress/",sep="") #define stressed model path
load_stress_files <- fs::dir_ls(path = stress_path, regexp = paste0(stress_path, "motif_.*\\.csv$"))

unstressed_path <- paste("Data/simulations/",n_spp,"_spp/",model,"/unstressed/",sep="") #define unstressed model path
load_unstressed_files <- fs::dir_ls(path = unstressed_path, regexp = paste0(unstressed_path, "motif_.*\\.csv$"))

raw_stress_data <- as.data.table(vroom::vroom(load_stress_files[grepl(paste0("_",motif,"\\."),load_stress_files)],show_col_types = FALSE)) %>%
  dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column")) #load stressed models and rename species columns

raw_unstressed_data <- as.data.table(vroom::vroom(load_unstressed_files[grepl(paste0("_",motif,"\\."),load_unstressed_files)],show_col_types = FALSE)) %>%
  dplyr::rename_with(~paste("spp",gsub("Column","",.x),sep = "_"),dplyr::starts_with("Column")) #load unstressed models and rename species columns

inflection_detection_stress <- julia_wrangle_dt(data = copy(raw_stress_data),
                                                n_spp = n_spp,
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
  merge(inflection_detection_stress,by="sim_id") %>%
  .[,first_neg := min(ifelse(apply(.SD, 1,min) < 0 , time, max(time))),.SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #due to erroneous negative values from simulation solver, must crop prior to first negative value
  .[.[,.I[time < first_neg], by = "sim_id"]$V1] %>% #filter prior to first negative value 
  .[,first_neg:=NULL] %>%
  .[sim_id %in% simulation_id,] %>%
  .[,stressed := "Stressed"]

resilience_unstressed_csv <- copy(raw_unstressed_data) %>% #repeat the inflection point analysis but for the unstressed models
  setnames(.,"target_spp","target_node") %>%
  .[,target_node := NA] %>% #as no node experiences stress, set to NA
  .[,community_id := paste(motif,community,sep = "_")] %>%
  .[,sim_id := paste(motif,community,sim,sep = "_")] %>%
  .[,time := seq_along(spp_1),by = "sim_id"] %>%
  .[,first_neg := min(ifelse(apply(.SD, 1,min) < 0 , time, max(time))),.SDcols = grep("^spp_",names(.)), by = "sim_id"] %>% #due to erroneous negative values from simulation solver, must crop prior to first negative value
  .[.[,.I[time < first_neg], by = "sim_id"]$V1] %>% #filter prior to first negative value 
  .[,first_neg:=NULL] %>% 
  .[,inflection_pt_incr := max(time), by = "sim_id"] %>% #as no stress, set inflection point to last time point for downstreaming processing
  .[,inflection_pt_decr :=  max(time), by = "sim_id"] %>%
  .[sim_id %in% simulation_id,] %>%
  .[,stressed := "Unstressed"]

rm(inflection_detection_stress,raw_stress_data,raw_unstressed_data) #remove object to release memory 

simulation_data <- rbind(resilience_stress_csv,resilience_unstressed_csv) %>%
  melt(.,measure.vars = names(.)[grep("^spp_",names(.))],
       variable.name = "species",value.name = "density")

#######################
# Create figure
#######################

p1 <- ggraph(example_graph, layout = "stress") +
  geom_edge_link(aes(color = shared),
                 width = 2) +
  geom_node_point(size = 11,
                  col = "#5D369D") +
  geom_edge_loop(aes(strength = 0.5),check_overlap = T, 
                 arrow = arrow(length = unit(6, "pt"), type = "closed"),
                 start_cap = circle(3, "mm"),
                 end_cap = circle(3, "mm"))  +
  scale_edge_color_manual(values = c("#E86100","#929292"),guide = "none")+
  geom_node_text(aes(label = species),size = 5, col="white") +
  scale_x_reverse() + 
  coord_flip(ylim = c(-0.75,0.8),xlim = c(2,-2)) +
  theme_void()


p2 <- ggplot(simulation_data,aes(x = time, y = density,col=species)) +
  geom_segment(aes(x = inflection_pt_incr,xend = inflection_pt_incr,y = 0, yend = Inf),
               col = "black",linetype="dashed") +
  geom_line() +
  facet_wrap(~stressed) +
  xlab("Time point") + ylab("Density") +
  scale_color_manual(values = c("grey40","#bfbd3d", "#5d3099", "#69c756","#6886c4"),
                     guide = "none") +
  theme_bw() +
  theme(legend.title = ggplot2::element_blank())

p3 <- ggplot() +
  ggpubr::background_image(magick::image_read("Results/figures/figure1_c.png"))

ggsave("Results/figures/figure1.pdf",
       (p1|p2)/p3  + 
         plot_annotation(tag_levels = "a"),
       width = 6,height = 6)
