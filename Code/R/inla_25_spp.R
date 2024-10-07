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

prior_prec = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"
ar_prior = list(theta1 = list(prior="pc.prec", param=c(0.01*3, 0.01)), #penalised complexity prior statement, where we state the odds of the standard deviation captured by the ar1 term exceeding 3*residual standard deviation are very small (probabaility = 0.01)
                theta2 = list(prior="pc.cor1", param=c(0, 0.9), initial = 0)) #Here we state the odds of detecting some autocorrelation (rho greater than 0) are very high (probability of 0.9). We set the inital rho at 0.

pred_data <- expand.grid("stand_time" = seq(0,1,0.01), #range of data to visualise trends over
                         "stressed" = paste(c(0,1)),
                         "ts_length" = c(0.1,0.4,0.7),
                         "search_effort" = c(0.1,0.5,1.0))

draws <- 5000 #how samples from posterior for coefficient estimates
newdata <- expand.grid("stand_time" = c(0,1), #only provide max and min stand_time to estimate slope
                       "stressed" = paste(c(0,1)),
                       "ts_length" = c(0.2,0.4,0.6),
                       "search_effort" = seq(0.1,1.0,0.2))
Xmat <- model.matrix(~stand_time*stressed*ts_length*search_effort, data=newdata)
lincomb <- inla.make.lincombs(as.data.frame(Xmat)) #linear combinations to provide to INLA

sample_id <- expand.grid("motif" = 1,
                         "community" = 1:25,
                         "sim" = 1:25) |>
  dplyr::mutate(sample_id = paste(motif,community,sim,sep = "_")) |>
  dplyr::select(sample_id) |> unlist() |> unname()

#prep_id <- readRDS("Data/resilience/sample_id.rds")

#######################
# Prepare 25 spp data
#######################

resilience_25_data <- as.data.table(readRDS("Data/resilience/full/motif1_25_invasive.rds")) %>%
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
# 25 spp multiJI
#######################

inla_data_25_1 <- subset(resilience_25_data,metric == "multiJI") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit()|>
  dplyr::select(metric_value,stand_time,stressed,ts_length,search_effort,comm_id,ar_id,ar_id_numeric,time)

#Zcomm_sim <- as(model.matrix( ~ 0 + comm_id:sim_id, data = inla_data_25_1), "Matrix") #unable to fit nested effect due to memory limitations

#fit model
inla_25_1 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                  +
                    f(comm_id, model = "iid",
                      hyper = prior_prec)  #random intercept for shared communities
                  +f(ar_id,stand_time, model = "iid",
                     hyper =  prior_prec)  #random slopes per simulation
                  +f(time, model = "ar",order = 1,replicate = ar_id_numeric,hyper = ar_prior) #autocorrelation term
                  ,
                  data = inla_data_25_1, 
                  family = "gaussian",
                  lincomb = lincomb,
                  control.inla = list(
                    int.strategy = "eb",
                    force.diagonal = T
                  ),
                  num.threads = 4,
                  #control.predictor=list(compute=TRUE),
                  #control.compute = list(config = TRUE,cpo = F),
                  control.fixed = prior.fixed
)

plot((1:nrow(na.omit(inla_data_25_1)))/(nrow(na.omit(inla_data_25_1)) +1), 
     sort(na.omit(inla_25_1$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1) #expect points to be approximately uniformly distributed across line
any(na.omit(inla_25_1$cpo$failure)>=1) #expect 0/few cpo >=1

marg_25_1 <- inla_25_1$marginals.lincomb.derived[grepl("lc",names(inla_25_1$marginals.lincomb.derived))] #subset marginals (equivalet to posteriors) to the linear combinations
groups_25_1 <- rep(1:(length(marg_25_1)/2),each = 2) #assign groups to match newdata (two posteriors per group; start and end of stand_time)
marg_grouped_25_1 <- split(marg_25_1,groups_25_1) #split the marginal list by these groups

slope_post_25_1 <- lapply(marg_grouped_25_1,function(x){
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

slope_range_25_1 <- slope_post_25_1 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5)) #summarise the slope posteriors

ggplot(slope_range_25_1,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_25_1, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  #coord_cartesian(ylim = c(-0.2,0.35))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_25_1,"Results/models/motif1_25_invasive_multiJI_model.rds")

#######################
# 25 spp mean_uniJI
#######################

inla_data_25_2 <-  subset(resilience_25_data,metric == "mean_uniJI" ) |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit() |>
  dplyr::select(metric_value,stand_time,stressed,ts_length,search_effort,comm_id,ar_id,ar_id_numeric,time)

inla_25_2 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                  +
                    f(comm_id, model = "iid",
                      hyper = prior_prec)  #random intercept for shared communities
                  +f(ar_id,stand_time, model = "iid",
                     hyper =  prior_prec)  #random slopes per simulation
                  +f(time, model = "ar",order = 1,replicate = ar_id_numeric,hyper = ar_prior) #autocorrelation term
                  ,
                  data = inla_data_25_2, 
                  family = "gaussian",
                  lincomb = lincomb,
                  control.inla = list(int.strategy = "eb",
                                      force.diagonal = T),
                  num.threads = 4,
                  #control.predictor=list(compute=TRUE),
                  #control.compute = list(config = TRUE,cpo = F),
                  control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_25_2)))/(nrow(na.omit(inla_data_25_2)) +1), 
     sort(na.omit(inla_25_2$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_25_2$cpo$failure)>=1)

marg_25_2 <- inla_25_2$marginals.lincomb.derived[grepl("lc",names(inla_25_2$marginals.lincomb.derived))]
groups_25_2 <- rep(1:(length(marg_25_2)/2),each = 2)
marg_grouped_25_2 <- split(marg_25_2,groups_25_2)

slope_post_25_2 <- lapply(marg_grouped_25_2,function(x){
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

slope_range_25_2 <- slope_post_25_2 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_25_2,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_25_2, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  coord_cartesian(ylim = c(-0.2,0.35))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_25_2,"Results/models/motif1_25_invasive_uniJI_model.rds")

#######################
# 25 spp FI
#######################

inla_data_25_3 <-  subset(resilience_25_data,metric == "FI") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit()|>
  dplyr::select(metric_value,stand_time,stressed,ts_length,search_effort,comm_id,ar_id,ar_id_numeric,time)

inla_25_3 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                  +
                    f(comm_id, model = "iid",
                      hyper = prior_prec)  #random intercept for shared communities
                  +f(ar_id,stand_time, model = "iid",
                     hyper =  prior_prec)  #random slopes per simulation
                  +f(time, model = "ar",order = 1,replicate = ar_id_numeric,hyper = ar_prior) #autocorrelation term
                  ,
                  data = inla_data_25_3, 
                  family = "gaussian",
                  lincomb = lincomb,
                  control.inla = list(int.strategy = "eb",
                                      force.diagonal = T),
                  num.threads = 4,
                  #control.predictor=list(compute=TRUE),
                  #control.compute = list(config = TRUE,cpo = F),
                  control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_25_3)))/(nrow(na.omit(inla_data_25_3)) +1), 
     sort(na.omit(inla_25_3$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_25_3$cpo$failure)>=1)

marg_25_3 <- inla_25_3$marginals.lincomb.derived[grepl("lc",names(inla_25_3$marginals.lincomb.derived))]
groups_25_3 <- rep(1:(length(marg_25_3)/2),each = 2)
marg_grouped_25_3 <- split(marg_25_3,groups_25_3)

slope_post_25_3 <- lapply(marg_grouped_25_3,function(x){
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

slope_range_25_3 <- slope_post_25_3 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_25_3,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_25_3, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  #coord_cartesian(ylim = c(-3.5,1.0))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_25_3,"Results/models/motif1_25_invasive_FI_model.rds")

#######################
# 25 spp MVI
#######################

inla_data_25_4 <-  subset(resilience_25_data,metric == "mvi") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  dplyr::mutate(metric_value = log(metric_value)) |>
  na.omit() |>
  dplyr::select(metric_value,stand_time,stressed,ts_length,search_effort,comm_id,ar_id,ar_id_numeric,time)

inla_25_4 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                  +
                    f(comm_id, model = "iid",
                      hyper = prior_prec)  #random intercept for shared communities
                  +f(ar_id,stand_time, model = "iid",
                     hyper =  prior_prec)  #random slopes per simulation
                  +f(time, model = "ar",order = 1,replicate = ar_id_numeric,hyper = ar_prior) #autocorrelation term
                  ,
                  data = inla_data_25_4, 
                  family = "gaussian",
                  lincomb = lincomb,
                  control.inla = list(int.strategy = "eb",
                                      force.diagonal = T),
                  num.threads = 4,
                  #control.predictor=list(compute=TRUE),
                  #control.compute = list(config = TRUE,cpo = F),
                  control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_25_4)))/(nrow(na.omit(inla_data_25_4)) +1), 
     sort(na.omit(inla_25_4$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_data_25_4$cpo$failure)>=1)

marg_25_4 <- inla_25_4$marginals.lincomb.derived[grepl("lc",names(inla_25_4$marginals.lincomb.derived))]
groups_25_4 <- rep(1:(length(marg_25_4)/2),each = 2)
marg_grouped_25_4 <- split(marg_25_4,groups_25_4)

slope_post_25_4 <- lapply(marg_grouped_25_4,function(x){
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

slope_range_25_4 <- slope_post_25_4 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_25_4,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_25_4, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  coord_cartesian(ylim = c(-0.1,0.7))+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_25_4,"Results/models/motif1_25_invasive_MVI_model.rds")

#######################
# 25 spp multiAR
#######################

inla_data_25_5 <-  subset(resilience_25_data,metric == "multiAR") |>
  dplyr::mutate(stressed = paste(stressed),
                ar_id_numeric = as.numeric(as.factor(ar_id))) |>
  na.omit() |>
  dplyr::select(metric_value,stand_time,stressed,ts_length,search_effort,comm_id,ar_id,ar_id_numeric,time)

inla_25_5 <- inla(metric_value ~ stand_time*stressed*ts_length*search_effort 
                  +
                    f(comm_id, model = "iid",
                      hyper = prior_prec)  #random intercept for shared communities
                  +f(ar_id,stand_time, model = "iid",
                     hyper =  prior_prec)  #random slopes per simulation
                  +f(time, model = "ar",order = 1,replicate = ar_id_numeric,hyper = ar_prior) #autocorrelation term
                  ,
                  data = inla_data_25_5, 
                  family = "gaussian",
                  lincomb = lincomb,
                  control.inla = list(int.strategy = "eb",
                                      force.diagonal = T),
                  num.threads = 4,
                  #control.predictor=list(compute=TRUE),
                  #control.compute = list(config = TRUE,cpo = F),
                  control.fixed = prior.fixed)

plot((1:nrow(na.omit(inla_data_25_5)))/(nrow(na.omit(inla_data_25_5)) +1), 
     sort(na.omit(inla_25_5$cpo$pit)),
     xlab="Uniform quantiles", ylab="Sorted PIT values")
abline(0,1)
any(na.omit(inla_data_25_5$cpo$failure)>=1)

marg_25_5 <- inla_25_5$marginals.lincomb.derived[grepl("lc",names(inla_25_5$marginals.lincomb.derived))]
groups_25_5 <- rep(1:(length(marg_25_5)/2),each = 2)
marg_grouped_25_5 <- split(marg_25_5,groups_25_5)

slope_post_25_5 <- lapply(marg_grouped_25_5,function(x){
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

slope_range_25_5 <- slope_post_25_5 |>
  dplyr::group_by(stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

ggplot(slope_range_25_5,aes(y = .value, x = search_effort,group = stressed)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= slope_post_25_5, aes(fill=stressed,group = stressed),alpha=0.5,position = position_dodge(width=0.1)) +
  tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=stressed,),
                                interval_size_range = c(0.8, 2),position = position_dodge(width=0.1)) +
  scale_color_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"),guide = "none")+
  scale_fill_manual(values = c("#67705F","#E3A59F"),name = "Presence\nof stress",labels = c("No","Yes"))+
  facet_wrap(~ts_length)+
  xlab("Search effort") +
  ylab("Trend estimate") + 
  theme_bw()

saveRDS(inla_25_5,"Results/models/motif1_25_invasive_multiAR_model.rds")

#######################
# Save 25 spp slopes
#######################
save(slope_post_25_1,slope_post_25_2,slope_post_25_3,slope_post_25_4,slope_post_25_5,
     file = "Results/models/motif1_25_invasive_slope_posteriors.RData")

save(slope_range_25_1,slope_range_25_2,slope_range_25_3,slope_range_25_4,slope_range_25_5,
     file = "Results/models/motif1_25_invasive_slope_ranges.RData")
