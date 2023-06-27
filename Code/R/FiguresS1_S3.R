require(data.table)
require(magrittr)
require(dplyr)
require(ggplot2)
require(patchwork)

#######################
# Load spp data
#######################
set.seed(123)

sample_id <- expand.grid("motif" = 1,
                         "community" = 1:30,
                         "sim" = 1:25) |>
  dplyr::mutate(sample_id = paste(motif,community,sim,sep = "_")) |>
  dplyr::select(sample_id) |> unlist() |> unname() |>
  sample(10)

resilience_5_data <- as.data.table(readRDS("Data/resilience/full/motif1_5_invasive.rds")) %>%
  .[sim_id %in% sample_id,] %>%
  .[,comm_id := paste0(unlist(strsplit(sim_id, "_", fixed=TRUE))[1:2],collapse = "_"), 
    by = sim_id] %>% #find shared community id
  .[,c("sim_id","comm_id","time","multiJI","mean_uniJI","FI","mvi","model","stressed","ts_length","search_effort","motif")] %>%
  #.[,mvi := log(mvi)] %>%
  .[,ar_id := paste(sim_id,model,stressed,ts_length,search_effort,sep = "_")] %>%
  melt(.,measure.vars = c("multiJI","mean_uniJI","FI","mvi"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot_loner
  .[search_effort == 1.0,] %>%
  .[,n_spp := 5] 

resilience_15_data <- as.data.table(readRDS("Data/resilience/full/motif1_15_invasive.rds")) %>%
  .[sim_id %in% sample_id,] %>%
  .[,comm_id := paste0(unlist(strsplit(sim_id, "_", fixed=TRUE))[1:2],collapse = "_"), 
    by = sim_id] %>% #find shared community id
  .[,c("sim_id","comm_id","time","multiJI","mean_uniJI","FI","mvi","model","stressed","ts_length","search_effort","motif")] %>%
  #.[,mvi := log(mvi)] %>%
  .[,ar_id := paste(sim_id,model,stressed,ts_length,search_effort,sep = "_")] %>%
  melt(.,measure.vars = c("multiJI","mean_uniJI","FI","mvi"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot_loner
  .[search_effort == 1.0,] %>%
  .[,n_spp := 15] 

resilience_25_data <- as.data.table(readRDS("Data/resilience/full/motif1_25_invasive.rds")) %>%
  .[sim_id %in% sample_id,] %>%
  .[,comm_id := paste0(unlist(strsplit(sim_id, "_", fixed=TRUE))[1:2],collapse = "_"), 
    by = sim_id] %>% #find shared community id
  .[,c("sim_id","comm_id","time","multiJI","mean_uniJI","FI","mvi","model","stressed","ts_length","search_effort","motif")] %>%
  #.[,mvi := log(mvi)] %>%
  .[,ar_id := paste(sim_id,model,stressed,ts_length,search_effort,sep = "_")] %>%
  melt(.,measure.vars = c("multiJI","mean_uniJI","FI","mvi"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot_loner
  .[search_effort == 1.0,] %>%
  .[,n_spp := 25] 

figure_data <- rbind(resilience_5_data,resilience_15_data,resilience_25_data) %>%
  .[,stand_time := scales::rescale(time),
    by = c("ar_id","metric","n_spp")] %>% #rescale time for comparability between 0-1
  .[,stressed := paste(stressed)] %>% #convert to factor
  .[,sim_id := paste0(unlist(strsplit(sim_id, "_", fixed=TRUE))[2:3],collapse = "_")] #trim shared motif suffix

#######################
# Create Figure S1
#######################
ggsave("Results/figures/figureS1.pdf",
ggplot(subset(figure_data,metric == "multiJI" & n_spp == 5),
       aes(x=stand_time, y= metric_value,col=stressed)) +
  geom_line() +
  facet_grid(ts_length~sim_id,scales = "free_y") +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
  scale_x_continuous(breaks = c(0.2,0.7))+
  xlab("Standardised time [0-1]") +
  ylab("multiJI metric value") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+

ggplot(subset(figure_data,metric == "mean_uniJI" & n_spp == 5),
       aes(x=stand_time, y= metric_value,col=stressed)) +
  geom_line() +
  facet_grid(ts_length~sim_id,scales = "free_y") +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
  scale_x_continuous(breaks = c(0.2,0.7))+
  xlab("Standardised time [0-1]") +
  ylab("mean_uniJI metric value") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank())+

ggplot(subset(figure_data,metric == "FI" & n_spp == 5),
       aes(x=stand_time, y= metric_value,col=stressed)) +
  geom_line() +
  facet_grid(ts_length~sim_id) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
  scale_x_continuous(breaks = c(0.2,0.7))+
  xlab("Standardised time [0-1]") +
  ylab("FI metric value") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+

ggplot(subset(figure_data,metric == "mvi" & n_spp == 5),
       aes(x=stand_time, y= metric_value,col=stressed)) +
  geom_line() +
  facet_grid(ts_length~sim_id,scales = "free_y") +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
  scale_x_continuous(breaks = c(0.2,0.7))+
  xlab("Standardised time [0-1]") +
  ylab("mvi metric value") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  plot_layout(nrow = 2,ncol=2, byrow = FALSE,guides = "collect"),
width = 12,height = 7)

#######################
# Create Figure S2
#######################
ggsave("Results/figures/figureS2.pdf",
       ggplot(subset(figure_data,metric == "multiJI" & n_spp == 15),
              aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id,scales = "free_y") +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("multiJI metric value") +
         theme_bw() +
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot(subset(figure_data,metric == "mean_uniJI" & n_spp == 15),
                aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id,scales = "free_y") +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("mean_uniJI metric value") +
         theme_bw() +
         theme(panel.grid.minor = element_blank(),
               panel.background = element_blank())+
         
         ggplot(subset(figure_data,metric == "FI" & n_spp == 15),
                aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id) +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("FI metric value") +
         theme_bw() +
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot(subset(figure_data,metric == "mvi" & n_spp == 15),
                aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id,scales = "free_y") +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("mvi metric value") +
         theme_bw() + 
         theme(panel.grid.minor = element_blank(),
               panel.background = element_blank()) +
         plot_layout(nrow = 2,ncol=2, byrow = FALSE,guides = "collect"),
       width = 12,height = 7)

#######################
# Create Figure S3
#######################
ggsave("Results/figures/figureS3.pdf",
       ggplot(subset(figure_data,metric == "multiJI" & n_spp == 25),
              aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id,scales = "free_y") +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("multiJI metric value") +
         theme_bw() +
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot(subset(figure_data,metric == "mean_uniJI" & n_spp == 25),
                aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id,scales = "free_y") +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("mean_uniJI metric value") +
         theme_bw() +
         theme(panel.grid.minor = element_blank(),
               panel.background = element_blank())+
         
         ggplot(subset(figure_data,metric == "FI" & n_spp == 25),
                aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id) +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("FI metric value") +
         theme_bw() +
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot(subset(figure_data,metric == "mvi" & n_spp == 25),
                aes(x=stand_time, y= metric_value,col=stressed)) +
         geom_line() +
         facet_grid(ts_length~sim_id,scales = "free_y") +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes")) +
         scale_x_continuous(breaks = c(0.2,0.7))+
         xlab("Standardised time [0-1]") +
         ylab("mvi metric value") +
         theme_bw() + 
         theme(panel.grid.minor = element_blank(),
               panel.background = element_blank()) +
         plot_layout(nrow = 2,ncol=2, byrow = FALSE,guides = "collect"),
       width = 12,height = 7)
