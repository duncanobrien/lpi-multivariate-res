# Graphical abstract

require(dplyr)
require(ggplot2)
require(patchwork)
require(tidybayes)
require(data.table)
require(magrittr)
source("Code/R/auxiliary_functions/julia_wrangle_function.R")
load("Results/models/motif1_5_invasive_slope_posteriors.RData")
load("Results/models/motif1_5_invasive_slope_ranges.RData")

set.seed(124)
motif = 1 
model = "invasive"
n_spp = 5
simulation_id <- "1_2_1" #example simulation

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
       variable.name = "species",value.name = "density") %>%
  .[,density := c(scale(density)),by = c("species","stressed")] %>%
  .[,density := density - density[1],by = c("species","stressed")]

#######################
# Create hypothesis figures
#######################
stress_df <- data.frame(x.poly = seq(100,200,0.1),
                        y.poly = seq(3.5,4.0, length.out = length(seq(100,200,0.1))),
                        xend=100,
                        yend=3.5) |>
  dplyr::mutate(stressed = "Stressed")

ts_df <- expand.grid(ts_length = c(20,40,60),
                     xmax = unique(simulation_data$inflection_pt_incr)) |>
  dplyr::mutate(stressed = ifelse(xmax >150,"Unstressed","Stressed"),
                xmin = xmax - ts_length,
                ymin = 4, 
                ymax = 4.5,
                id = row_number()) |>
  dplyr::mutate(labx = xmin + 10,
                laby = ymax - ((ymax-ymin)/2), .by = c(ts_length,stressed))

ggplot2::ggsave("/Users/ul20791/Desktop/Academia/Papers/OBrien_et_al_resilience_metrics/g_abstract_fig1.pdf",
                ggplot(simulation_data) +
                  geom_segment(aes(x = inflection_pt_incr,xend = inflection_pt_incr,y = -Inf, yend = Inf),
                               col = "black",linetype="dashed") +
                  geom_line(aes(x = time, y = density,col=species)) +
                  scale_color_manual(values = c("grey40","#bfbd3d", "#5d3099", "#69c756","#6886c4"),
                                     guide = "none") +
                  ggnewscale::new_scale_colour()+
                  geom_rect(data=ts_df, 
                            mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, group = id, fill = factor(ts_length)), alpha=0.4) +
                  geom_text(data=ts_df, 
                            mapping=aes(x=labx, y= laby, label = ts_length, group = id), col = "white") +
                  #scale_fill_manual(values = c("#EA4C46","#F07470", "#F1959B"),name = "Time series\nlength") +
                  scale_fill_manual(values = c("grey90","grey60", "grey20"),name = "Time series\nlength") +
                  facet_wrap(~stressed) +
                  xlab("Time point") + ylab("Scaled population size") +
                  ggnewscale::new_scale_colour()+
                  geom_segment(data = stress_df,
                               aes(x = x.poly, y = yend, xend = x.poly, yend = y.poly, color = y.poly))+
                  scale_colour_gradientn(colours = c("grey", "grey20"),guide = 'none') +
                  # annotate("rect", fill = "red", alpha = 0.5, 
                  #          xmin = ts_df$xmin, xmax = ts_df$xmax,
                  #          ymin = -Inf, ymax = Inf) +
                  theme_bw() +
                  theme(panel.grid= element_blank(),
                        strip.background = element_rect(fill="white"))
                ,
                width = 7, height = 3
)


ggplot2::ggsave("/Users/ul20791/Desktop/Academia/Papers/OBrien_et_al_resilience_metrics/g_abstract_fig2.pdf",
                ggplot2::ggplot(slope_range_5_1 |> mutate(metric = "multiJI") |>
                  mutate(ts_length = paste(ts_length,"years"),
                         search_effort = search_effort*100),
                aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="black") +
  geom_line(aes(col=stressed,group = stressed)) +
  tidybayes::stat_slab(data= slope_post_5_1 |> mutate(metric = "multiJI") |>
                         mutate(ts_length = paste(ts_length,"years"),
                                search_effort = search_effort*100),
                       aes(fill=stressed,group = stressed),
                       linewidth = 1,
                       alpha = 0.3,
                       position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed),
                                fatten_point = 1.3,
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#047101","#FF8EC6"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#047101","#FF8EC6"),name = "Presence\nof stress",labels = c("No","Yes"))+
  scale_alpha_manual(values = c(0.1,0.3),guide = "none") + 
  scale_x_continuous(breaks = seq(0.1,1.0,0.2)*100)+
  facet_wrap(~ts_length)+
  coord_cartesian(ylim = c(-0.1,0.35))+
  xlab("Search effort (%)") +
  ylab("Trend estimate") + 
  theme_bw()+
  theme(panel.grid= element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill="white"),
        plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt")),
  width = 5, height = 3
)


