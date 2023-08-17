require(data.table)
require(magrittr)
require(dplyr)
require(ggplot2)
require(INLA)

#######################
# Preamble
#######################
prior.fixed <- list(mean.intercept = 0, prec.intercept = 1,
                    mean = 0, prec = 1,
                    correlation.matrix = F) #normal prior on intercept and fixed effects (mean = 0, sd = 1)

pred_data <- expand.grid("stand_time" = seq(0,1,0.01), #range of data to visualise trends over
                         "stressed" = paste(c(0,1)),
                         "ts_length" = c(0.1,0.4,0.7),
                         "search_effort" = c(0.1,0.5,1.0))

draws <- 10000 #how samples from posterior for coefficient estimates
newdata <- expand.grid("stand_time" = c(0,1), #only provide max and min stand_time to estimate slope
                       "stressed" = paste(c(0,1)),
                       "ts_length" = c(0.2,0.4,0.6),
                       "search_effort" = seq(0.1,1.0,0.2))
Xmat <- model.matrix(~stand_time*stressed*ts_length*search_effort, data=newdata)
lincomb <- inla.make.lincombs(as.data.frame(Xmat)) #linear combinations to provide to INLA

sample_id <- expand.grid("motif" = 1,
                         "community" = 1:30,
                         "sim" = 1:25) |>
  dplyr::mutate(sample_id = paste(motif,community,sim,sep = "_")) |>
  dplyr::select(sample_id) |> unlist() |> unname()

#prep_id <- readRDS("Data/resilience/sample_id.rds")

#######################
# Prepare 15 spp data
#######################

resilience_15_data <- as.data.table(readRDS("Data/resilience/full/motif1_15_invasive.rds")) %>%
  .[sim_id %in% sample_id,] %>%
  .[,comm_id := paste0(unlist(strsplit(sim_id, "_", fixed=TRUE))[1:2],collapse = "_"), 
    by = sim_id] %>% #find shared community id
  .[,c("sim_id","comm_id","time","multiJI","mean_uniJI","FI","mvi","multiAR","model","stressed","ts_length","search_effort","motif")] %>%
  #.[,mvi := log(mvi)] %>%
  .[,ar_id := paste(sim_id,model,stressed,ts_length,search_effort,sep = "_")] %>%
  melt(.,measure.vars = c("multiJI","mean_uniJI","FI","mvi","multiAR"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot_loner
  .[,ts_length := ts_length/100] %>% #so on same scale as search_effort
  .[,stand_time := scales::rescale(time),
    by = c("ar_id","metric")] %>% #to prevent model extapolating ts_len 10 to length 70, ensure time is scaled
  # .[,metric_value := scale(metric_value),
  #   by = c("ts_length","search_effort","stressed","id","metric")] %>% #mvi and FI on different scale to JIs
  as.data.frame()

#######################
# 15 spp multiJI
#######################

inla_data_15_1 <- subset(resilience_15_data,metric == "multiJI") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit()

#fit model
inla_15_1 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                 +
                   f(comm_id, model = "iid",
                     hyper = list(prec = list(prior = "logtnormal",
                                              param = c(0, 1))))  #random intercept for shared communities
                 +f(ar_id,stand_time, model = "iid",
                    hyper =list(prec = list(prior = "logtnormal",
                                            param = c(0, 1))))  #random slopes per simulation
                 +f(time, model = "ar",order = 1,replicate = ar_id_numeric) #autocorrelation term
                 ,
                 data = inla_data_15_1, 
                 family = "gaussian",
                 lincomb = lincomb,
                 control.inla = list(
                   int.strategy = "eb",
                   force.diagonal = T
                 ),
                 num.threads = 5,
                 control.predictor=list(compute=TRUE),
                 control.compute = list(config = TRUE,cpo = F),
                 control.fixed = prior.fixed
)

plot((1:nrow(na.omit(inla_data_15_1)))/(nrow(na.omit(inla_data_15_1)) +1), 
     sort(na.omit(inla_15_1$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1) #expect points to be approximately uniformly distributed across line
any(na.omit(inla_15_1$cpo$failure)>=1) #expect 0/few cpo >=1

marg_15_1 <- inla_15_1$marginals.lincomb.derived[grepl("lc",names(inla_15_1$marginals.lincomb.derived))] #subset marginals (equivalet to posteriors) to the linear combinations
groups_15_1 <- rep(1:(length(marg_15_1)/2),each = 2) #assign groups to match newdata (two posteriors per group; start and end of stand_time)
marg_grouped_15_1 <- split(marg_15_1,groups_15_1) #split the marginal list by these groups

slope_post_15_1 <- lapply(marg_grouped_15_1,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x[[2]]) - 
               inla.rmarginal(draws, x[[1]])) #extract slope as difference between end posterior and start posterior
}) |>
  data.table::rbindlist() |>
  cbind(
    subset(newdata,stand_time == 1)[,c("stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws)) #merge with newdata
  ) |>
  dplyr::mutate(ts_length = ts_length*100) #convert back to true ts length for visualisations

slope_range_15_1 <- slope_post_15_1 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5)) #summarise the slope posteriors

ggplot(slope_range_15_1,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_15_1, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  #coord_cartesian(ylim = c(-0.2,0.35))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_15_1,"Results/models/motif1_15_invasive_multiJI_model.rds")

#######################
# 15 spp mean_uniJI
#######################

inla_data_15_2 <-  subset(resilience_15_data,metric == "mean_uniJI" ) |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit()

inla_15_2 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                 +
                   f(comm_id, model = "iid",
                     hyper = list(prec = list(prior = "logtnormal",
                                              param = c(0, 1)))) +
                   f(ar_id,stand_time, model = "iid",
                     hyper = list(prec = list(prior = "logtnormal",
                                              param = c(0, 1)))) 
                 +f(time, model = "ar",order = 1,replicate = ar_id_numeric) #autocorrelation term
                 ,
                 data = inla_data_15_2, 
                 family = "gaussian",
                 lincomb = lincomb,
                 control.inla = list(int.strategy = "eb",
                                     force.diagonal = T),
                 num.threads = 5,
                 control.predictor=list(compute=TRUE),
                 control.compute = list(config = TRUE,cpo = F),
                 control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_15_2)))/(nrow(na.omit(inla_data_15_2)) +1), 
     sort(na.omit(inla_15_2$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_15_1$cpo$failure)>=1)

marg_15_2 <- inla_15_2$marginals.lincomb.derived[grepl("lc",names(inla_15_2$marginals.lincomb.derived))]
groups_15_2 <- rep(1:(length(marg_15_2)/2),each = 2)
marg_grouped_15_2 <- split(marg_15_2,groups_15_2)

slope_post_15_2 <- lapply(marg_grouped_15_2,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x[[2]]) - 
               inla.rmarginal(draws, x[[1]]))
}) |>
  data.table::rbindlist() |>
  cbind(
    subset(newdata,stand_time == 1)[,c("stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws))
  ) |>
  dplyr::mutate(ts_length = ts_length*100)

slope_range_15_2 <- slope_post_15_2 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_15_2,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_15_2, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  coord_cartesian(ylim = c(-0.2,0.35))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_15_2,"Results/models/motif1_15_invasive_uniJI_model.rds")

#######################
# 15 spp FI
#######################

inla_data_15_3 <-  subset(resilience_15_data,metric == "FI") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit()

inla_15_3 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                 +
                   f(comm_id, model = "iid",
                     hyper =  list(prec = list(prior = "logtnormal",
                                               param = c(0, 1)))) +
                   f(ar_id,stand_time, model = "iid",
                     hyper =  list(prec = list(prior = "logtnormal",
                                               param = c(0, 1)))) 
                 +f(time, model = "ar",order = 1,replicate = ar_id_numeric) #autocorrelation term
                 ,
                 data = inla_data_15_3, 
                 family = "gaussian",
                 lincomb = lincomb,
                 control.inla = list(int.strategy = "eb",
                                     force.diagonal = T),
                 num.threads = 5,
                 control.predictor=list(compute=TRUE),
                 control.compute = list(config = TRUE,cpo = F),
                 control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_15_3)))/(nrow(na.omit(inla_data_15_3)) +1), 
     sort(na.omit(inla_15_3$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_15_3$cpo$failure)>=1)

marg_15_3 <- inla_15_3$marginals.lincomb.derived[grepl("lc",names(inla_15_3$marginals.lincomb.derived))]
groups_15_3 <- rep(1:(length(marg_15_3)/2),each = 2)
marg_grouped_15_3 <- split(marg_15_3,groups_15_3)

slope_post_15_3 <- lapply(marg_grouped_15_3,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x[[2]]) - 
               inla.rmarginal(draws, x[[1]]))
}) |>
  data.table::rbindlist() |>
  cbind(
    subset(newdata,stand_time == 1)[,c("stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws))
  ) |>
  dplyr::mutate(ts_length = ts_length*100)

slope_range_15_3 <- slope_post_15_3 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_15_3,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_15_3, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  #coord_cartesian(ylim = c(-3.5,1.0))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_15_3,"Results/models/motif1_15_invasive_FI_model.rds")

#######################
# 15 spp MVI
#######################

inla_data_15_4 <-  subset(resilience_15_data,metric == "mvi") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  dplyr::mutate(metric_value = log(metric_value)) |>
  na.omit()

inla_15_4 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                 +
                   f(comm_id, model = "iid",
                     hyper = list(prec = list(prior = "logtnormal",
                                              param = c(0, 1)))) +
                   f(ar_id,stand_time, model = "iid",
                     hyper = list(prec = list(prior = "logtnormal",
                                              param = c(0, 1)))) 
                 +f(time, model = "ar",order = 1,replicate = ar_id_numeric) #autocorrelation term
                 ,
                 data = inla_data_15_4, 
                 family = "gaussian",
                 lincomb = lincomb,
                 control.inla = list(int.strategy = "eb",
                                     force.diagonal = T),
                 num.threads = 5,
                 control.predictor=list(compute=TRUE),
                 control.compute = list(config = TRUE,cpo = F),
                 control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_15_4)))/(nrow(na.omit(inla_data_15_4)) +1), 
     sort(na.omit(inla_15_4$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_data_15_4$cpo$failure)>=1)

marg_15_4 <- inla_15_4$marginals.lincomb.derived[grepl("lc",names(inla_15_4$marginals.lincomb.derived))]
groups_15_4 <- rep(1:(length(marg_15_4)/2),each = 2)
marg_grouped_15_4 <- split(marg_15_4,groups_15_4)

slope_post_15_4 <- lapply(marg_grouped_15_4,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x[[2]]) - 
               inla.rmarginal(draws, x[[1]]))
}) |>
  data.table::rbindlist() |>
  cbind(
    subset(newdata,stand_time == 1)[,c("stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws))
  ) |>
  dplyr::mutate(ts_length = ts_length*100)

slope_range_15_4 <- slope_post_15_4 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_15_4,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_15_4, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  coord_cartesian(ylim = c(-0.1,0.7))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_15_4,"Results/models/motif1_15_invasive_MVI_model.rds")

#######################
# 15 spp multiAR
#######################

inla_data_15_5 <-  subset(resilience_15_data,metric == "multiAR") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit()

inla_15_5 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                  +
                    f(comm_id, model = "iid",
                      hyper = list(prec = list(prior = "logtnormal",
                                               param = c(0, 1)))) +
                    f(ar_id,stand_time, model = "iid",
                      hyper = list(prec = list(prior = "logtnormal",
                                               param = c(0, 1)))) 
                  +f(time, model = "ar",order = 1,replicate = ar_id_numeric) #autocorrelation term
                  ,
                  data = inla_data_15_5, 
                  family = "gaussian",
                  lincomb = lincomb,
                  control.inla = list(int.strategy = "eb",
                                      force.diagonal = T),
                  num.threads = 5,
                  control.predictor=list(compute=TRUE),
                  control.compute = list(config = TRUE,cpo = F),
                  control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_15_5)))/(nrow(na.omit(inla_data_15_5)) +1), 
     sort(na.omit(inla_15_5$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_data_15_5$cpo$failure)>=1)

marg_15_5 <- inla_15_5$marginals.lincomb.derived[grepl("lc",names(inla_15_5$marginals.lincomb.derived))]
groups_15_5 <- rep(1:(length(marg_15_5)/2),each = 2)
marg_grouped_15_5 <- split(marg_15_5,groups_15_5)

slope_post_15_5 <- lapply(marg_grouped_15_5,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x[[2]]) - 
               inla.rmarginal(draws, x[[1]]))
}) |>
  data.table::rbindlist() |>
  cbind(
    subset(newdata,stand_time == 1)[,c("stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws))
  ) |>
  dplyr::mutate(ts_length = ts_length*100)

slope_range_15_5 <- slope_post_15_5 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_15_5,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_15_5, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  coord_cartesian(ylim = c(-0.1,0.7))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_15_5,"Results/models/motif1_15_invasive_multiAR_model.rds")

#######################
# Save 15 spp slopes
#######################
save(slope_post_15_1,slope_post_15_2,slope_post_15_3,slope_post_15_4,slope_post_15_5,
     file = "Results/models/motif1_15_invasive_slope_posteriors2.RData")

save(slope_range_15_1,slope_range_15_2,slope_range_15_3,slope_range_15_4,slope_range_15_5,
     file = "Results/models/motif1_15_invasive_slope_ranges2.RData")
