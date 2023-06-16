# Figure 2

require(dplyr)
require(ggplot2)
require(patchwork)
require(tidybayes)

load("Results/models/motif1_15_invasive_slope_posteriors.RData")
load("Results/models/motif1_15_invasive_slope_ranges.RData")

#######################
# Create figure
#######################
ggplot2::ggsave("Results/figures/figure2.pdf",
       
       ggplot(slope_range_15_1 |> mutate(metric = "multiJI"),
              aes(y = .value, x = search_effort,group = stressed)) +
         geom_hline(yintercept = 0, linetype = "dashed", colour="black") +
         tidybayes::stat_slab(data= slope_post_15_1 |> mutate(metric = "multiJI"),
                              aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
         tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed),
                                       fatten_point = 1.3,
                                       interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
         scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
         scale_x_continuous(breaks = seq(0.1,1.0,0.2))+
         facet_grid(metric~ts_length)+
         coord_cartesian(ylim = c(-0.15,0.35))+
         xlab("Search effort") +
         ylab("Trend estimate") + 
         theme_bw()+
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot( slope_range_15_2 |> mutate(metric = "mean_uniJI"),
                 aes(y = .value, x = search_effort,group = stressed)) +
         geom_hline(yintercept = 0, linetype = "dashed", colour="black") +
         tidybayes::stat_slab(data = slope_post_15_2 |> mutate(metric = "mean_uniJI"),
                              aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
         tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed),
                                       fatten_point = 1.3,
                                       interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
         scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
         scale_x_continuous(breaks = seq(0.1,1.0,0.2))+
         facet_grid(metric~ts_length)+
         coord_cartesian(ylim = c(-0.075,0.175))+
         xlab("Search effort") +
         ylab("Trend estimate") + 
         theme_bw()+
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot(slope_range_15_3 |> mutate(metric = "FI"),
                aes(y = .value, x = search_effort,group = stressed)) +
         geom_hline(yintercept = 0, linetype = "dashed", colour="black") +
         tidybayes::stat_slab(data = slope_post_15_3 |> mutate(metric = "FI"),
                              aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
         tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed),
                                       fatten_point = 1.3,
                                       interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
         scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
         scale_x_continuous(breaks = seq(0.1,1.0,0.2))+
         facet_grid(metric~ts_length)+
         coord_cartesian(ylim = c(-1.05,0.075))+
         xlab("Search effort") +
         ylab("Trend estimate") + 
         theme_bw() +
         theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               plot.margin = margin(5.5,5.5,0.5,5.5, unit = "pt"))+
         
         ggplot( slope_range_15_4 |> mutate(metric = "mvi"),
                 aes(y = .value, x = search_effort,group = stressed)) +
         geom_hline(yintercept = 0, linetype = "dashed", colour="black") +
         tidybayes::stat_slab(data =  slope_post_15_4 |> mutate(metric = "mvi"),
                              aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
         tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed),
                                       fatten_point = 1.3,
                                       interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
         scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
         scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
         scale_x_continuous(breaks = seq(0.1,1.0,0.2))+
         facet_grid(metric~ts_length)+
         coord_cartesian(ylim = c(-0.075,0.6))+
         xlab("Search effort") +
         ylab("Log trend estimate") + 
         theme_bw()+ 
         theme(panel.grid.minor = element_blank(),
               panel.background = element_blank())+
         plot_layout(nrow = 4, byrow = FALSE,guides = "collect"),
       width = 6,height = 7)