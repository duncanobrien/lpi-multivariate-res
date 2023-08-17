# Figure 4

require(dplyr)
require(ggplot2)
require(patchwork)
require(tidybayes)

load("Results/models/motif1_threshold_posteriors.RData")
load("Results/models/motif1_threshold_ranges.RData")

summary_range_data <- rbind(summary_range1 |> dplyr::mutate(metric = "multiJI"),
                            summary_range2 |> dplyr::mutate(metric = "max_uniJI"),
                            summary_range3 |> dplyr::mutate(metric = "multiAR")) |>
  dplyr::mutate(n_spp = factor(n_spp,levels = c(5,15,25)))

summary_post_data <- rbind(summary_post1 |> dplyr::mutate(metric = "multiJI"),
                           summary_post2 |> dplyr::mutate(metric = "max_uniJI"),
                           summary_post3 |> dplyr::mutate(metric = "multiAR")) |>
  dplyr::mutate(n_spp = factor(n_spp,levels = c(5,15,25)))

inv_logit_perc <- scales::trans_new("inv_logit_perc",
                                    transform = function(x){suppressWarnings(plogis(x))},
                                    inverse = function(x){suppressWarnings(qlogis(x))})
labels <- c(-4,round(qlogis(c(0.25,0.5,0.75)),1),4) #ensure breaks @ 0.25 probability intervals
breaks <- c(-4,round(qlogis(c(0.25,0.5,0.75)),1),4)

#######################
# Create figure
#######################

# ggplot2::ggsave("Results/figures/figure4.pdf",
# ggplot(subset(summary_range_data,n_spp == 5),aes(x=.value,y=metric,group = stressed)) +
#   geom_vline(xintercept = 0,color = "black",linetype = "dashed")+
#   tidybayes::stat_slab(data= subset(summary_post_data,n_spp == 5), aes(fill=stressed,group = stressed),alpha=0.5,
#                        position = ggstance::position_dodgev(height = 0.75)) +
#   tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper,col=stressed),
#                                 fatten_point = 1.1,
#                                 interval_size_range = c(0.7, 1.5),position = ggstance::position_dodgev(height = 0.75)) +
#   facet_grid(search_effort~ts_length,scales = "free_y") + 
#   scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress")+
#   scale_color_manual(values = c("#67705F","#E3A59F"),
#                      guide = "none")+
#   xlab("Probability of exceeding 1") + ylab("Metric") +
#   scale_x_continuous(labels = round(plogis(labels),1), breaks = breaks)+
#   scale_y_discrete(expand = expansion(add = c(0.75, 1)))+
#   coord_trans(x = inv_logit_perc,xlim = c(-4.5,4.5))+
#   ggtitle("N = 5") +
#   theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.title.y=element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.key.width = unit(1.5, "line"),
#         plot.margin = margin(5.5,5.5,0.0,5.5, unit = "pt")) +
#   
#   plot_spacer() +
#   
#   ggplot(subset(summary_range_data,n_spp == 15),aes(x=.value,y=metric,group = stressed)) +
#   geom_vline(xintercept = 0,color = "black",linetype = "dashed")+
#   tidybayes::stat_slab(data= subset(summary_post_data,n_spp == 15), aes(fill=stressed,group = stressed),alpha=0.5,
#                        position = ggstance::position_dodgev(height = 0.75)) +
#   tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper,col=stressed),
#                                 fatten_point = 1.1,
#                                 interval_size_range = c(0.7, 1.5),position = ggstance::position_dodgev(height = 0.75)) +
#   facet_grid(search_effort~ts_length,scales = "free_y") + 
#   scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress")+
#   scale_color_manual(values = c("#67705F","#E3A59F"),
#                      guide = "none")+
#   xlab("Probability of exceeding 1") + ylab("Metric") +
#   scale_x_continuous(labels = round(brms::inv_logit_scaled(labels),1), breaks = breaks)+
#   scale_y_discrete(expand = expansion(add = c(0.75, 1)))+
#   coord_trans(x = inv_logit_perc,xlim = c(-4.5,4.5))+
#   ggtitle("N = 15") +
#   theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.key.width = unit(1.5, "line"),
#         plot.margin = margin(5.5,5.5,0.0,5.5, unit = "pt")) +
#   
#   plot_spacer() +
#   
#   ggplot(subset(summary_range_data,n_spp == 25),aes(x=.value,y=metric,group = stressed)) +
#   geom_vline(xintercept = 0,color = "black",linetype = "dashed")+
#   tidybayes::stat_slab(data= subset(summary_post_data,n_spp == 25), aes(fill=stressed,group = stressed),alpha=0.5,
#                        position = ggstance::position_dodgev(height = 0.75)) +
#   tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper,col=stressed),
#                                 fatten_point = 1.1,
#                                 interval_size_range = c(0.7, 1.5),position = ggstance::position_dodgev(height = 0.75)) +
#   facet_grid(search_effort~ts_length,scales = "free_y") + 
#   scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress")+
#   scale_color_manual(values = c("#67705F","#E3A59F"),
#                      guide = "none")+
#   xlab("Probability of exceeding 1") + ylab("Metric") +
#   scale_x_continuous(labels = round(brms::inv_logit_scaled(labels),1), breaks = breaks)+
#   scale_y_discrete(expand = expansion(add = c(0.75, 1)))+
#   coord_trans(x = inv_logit_perc,xlim = c(-4.5,4.5))+
#   ggtitle("N = 25") +
#   theme_bw() +
#   theme(axis.title.y=element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.key.width = unit(1.5, "line")) +
#   
#   plot_annotation(tag_levels ="a") +
#   plot_layout(nrow = 5, guides = "collect", heights = c(1,-0.105,1,-0.105,1)) &
#   theme(plot.title = element_text(size = 10, face = "bold")),
# width = 7,height = 7)
# 


ggplot2::ggsave("Results/figures/figure4_alt.pdf",
                ggplot(subset(summary_range_data,metric == "multiJI"),aes(x=.value,y=n_spp,group = stressed)) +
                  geom_vline(xintercept = 0,color = "black",linetype = "dashed")+
                  tidybayes::stat_slab(data= subset(summary_post_data,metric == "multiJI"), 
                                       aes(fill=stressed),alpha=0.5,
                                       position = ggstance::position_dodgev(height = 0.5)) +
                  tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper,col=stressed),
                                                fatten_point = 1.1,
                                                interval_size_range = c(0.7, 1.5),position = ggstance::position_dodgev(height = 0.5)) +
                  facet_grid(search_effort~ts_length,scales = "free_y") + 
                  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress")+
                  scale_color_manual(values = c("#67705F","#E3A59F"),
                                     guide = "none")+
                  xlab("Probability of exceeding 1") + ylab("Community size (#species)") +
                  scale_x_continuous(labels = round(plogis(labels),1), breaks = breaks)+
                  coord_trans(x = inv_logit_perc,xlim = c(-4.5,4.5))+
                  ggtitle("multiJI") +
                  theme_bw() +
                  theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.key.width = unit(1.5, "line"),
                        plot.margin = margin(5.5,5.5,0.0,5.5, unit = "pt")) +
                  
                  ggplot(subset(summary_range_data,metric == "max_uniJI"),aes(x=.value,y=n_spp,group = stressed)) +
                  geom_vline(xintercept = 0,color = "black",linetype = "dashed")+
                  tidybayes::stat_slab(data= subset(summary_post_data,metric == "max_uniJI"), aes(fill=stressed,group = stressed),alpha=0.5,
                                       position = ggstance::position_dodgev(height = 0.5)) +
                  tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper,col=stressed),
                                                fatten_point = 1.1,
                                                interval_size_range = c(0.7, 1.5),position = ggstance::position_dodgev(height = 0.5)) +
                  facet_grid(search_effort~ts_length,scales = "free_y") + 
                  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress")+
                  scale_color_manual(values = c("#67705F","#E3A59F"),
                                     guide = "none")+
                  xlab("Probability of exceeding 1") + ylab("Community size (#species)") +
                  scale_x_continuous(labels = round(plogis(labels),1), breaks = breaks)+
                  coord_trans(x = inv_logit_perc,xlim = c(-4.5,4.5))+
                  ggtitle("max_uniJI") +
                  theme_bw() +
                  theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.key.width = unit(1.5, "line"),
                        plot.margin = margin(5.5,5.5,0.0,5.5, unit = "pt")) +
                  
                  ggplot(subset(summary_range_data,metric == "multiAR"),aes(x=.value,y=n_spp,group = stressed)) +
                  geom_vline(xintercept = 0,color = "black",linetype = "dashed")+
                  tidybayes::stat_slab(data= subset(summary_post_data,metric == "multiAR"), aes(fill=stressed,group = stressed),alpha=0.5,
                                       position = ggstance::position_dodgev(height = 0.5)) +
                  tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper,col=stressed),
                                                fatten_point = 1.1,
                                                interval_size_range = c(0.7, 1.5),position = ggstance::position_dodgev(height = 0.5)) +
                  facet_grid(search_effort~ts_length,scales = "free_y") + 
                  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress")+
                  scale_color_manual(values = c("#67705F","#E3A59F"),
                                     guide = "none")+
                  xlab("Probability of exceeding threshold") + ylab("Community size (#species)") +
                  scale_x_continuous(labels = round(plogis(labels),1), breaks = breaks)+
                  coord_trans(x = inv_logit_perc,xlim = c(-4.5,4.5))+
                  ggtitle("multiAR") +
                  theme_bw() +
                  theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.key.width = unit(1.5, "line")) +
                  
                  plot_annotation(tag_levels ="a") +
                  plot_layout(nrow = 3, byrow = F,guides = "collect") &
                  theme(plot.title = element_text(size = 10, face = "bold")),
                width = 7,height = 7)
