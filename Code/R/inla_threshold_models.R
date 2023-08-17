require(data.table)
require(magrittr)
require(dplyr)
require(ggplot2)
require(INLA)
require(brms)
require(patchwork)

#######################
# Preamble
#######################
prior.fixed <- list(mean.intercept = 0, prec.intercept = 1,
                    mean = 0, prec = 1)

pred_success_data <- expand.grid("n_spp" = c(5,15,25),
                                 "stressed" = paste(c(0,1)),
                                 "ts_length" = c(0.1,0.4,0.7),
                                 "search_effort" = c(0.1,0.5,1.0))

sample_id <- expand.grid("motif" = 1,
                         "community" = 1:30,
                         "sim" = 1:25) |>
  dplyr::mutate(sample_id = paste(motif,community,sim,sep = "_")) |>
  dplyr::select(sample_id) |> unlist() |> unname()

draws <- 10000 #how samples from posterior for coefficient estimates
newdata <- expand.grid("n_spp" = c(5,15,25),
                       "stressed" = paste(c(0,1)),
                       "ts_length" = c(0.2,0.4,0.6),
                       "search_effort" = c(0.25,0.5,0.75,1.0)
)
Xmat_summary <- model.matrix(~n_spp*stressed*ts_length*search_effort, data=newdata)
lincomb_summary <- inla.make.lincombs(as.data.frame(Xmat_summary)) #linear combinations to provide to INLA

#######################
# Prepare data
#######################
files <- paste0("Data/resilience/summary/summary_motif1_",c(5,15,25),"_invasive.csv")

summary_data <- vroom::vroom(files,col_select = !starts_with("uniJI")) %>%
  as.data.table() %>%
  .[sim_id %in% sample_id,] %>%
  .[metric %in% c("multiJI","max_uniJI","mean_uniJI","multiAR"),] %>%
  .[,comm_id := paste(n_spp,
                       paste0(unlist(strsplit(sim_id, "_", fixed=TRUE))[1:2],collapse = "_"),
                       sep="_"), 
    by = sim_id] %>% #find shared community id
  .[,c("sim_id","metric","comm_id","threshold_crossed","corrected_trend","model","stressed","ts_length","search_effort","n_spp","motif")] %>%
  .[,ts_length := ts_length/100] %>% #so on same scale as search_effort
  .[,stressed := paste(stressed)] %>%
  #.[,n_spp := factor(n_spp,levels = c("5","15","25"))] %>%
  as.data.frame()

#######################
# Fit multijI
#######################

inla_success_data1 <-  subset(summary_data,metric == "multiJI") |>
  dplyr::bind_rows(pred_success_data)

inla_success1 <- inla(threshold_crossed ~ n_spp*stressed*ts_length*search_effort 
                      +
                        f(comm_id, model = "iid",
                          hyper = list(prec = list(prior = "logtnormal",
                                                   param = c(0, 1)))) ,
                      data = inla_success_data1, 
                      family = "binomial",
                      lincomb = lincomb_summary,
                      control.compute = list(config = TRUE,cpo = TRUE),
                      control.predictor=list(compute=TRUE),
                      control.fixed = prior.fixed)

any(na.omit(inla_success1$cpo$failure)>=1)

marg_summary1 <- inla_success1$marginals.lincomb.derived[grepl("lc",names(inla_success1$marginals.lincomb.derived))] #subset marginals (equivalet to posteriors) to the linear combinations

summary_post1 <- lapply(marg_summary1,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x))
}) |>
  data.table::rbindlist() |>
  cbind(
    newdata[,c("n_spp","stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws)) #merge with newdata
  ) |>
  dplyr::mutate(ts_length = ts_length*100,
                stressed = ifelse(stressed == 1,"Yes","No")) #convert back to true ts length for visualisations

summary_range1 <- summary_post1 |>
  dplyr::group_by(n_spp,stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

#######################
# Fit max_uniJI
#######################

inla_success_data2 <- subset(summary_data,metric == "max_uniJI") %>%
  dplyr::bind_rows(pred_success_data)

inla_success2 <- inla(threshold_crossed ~ n_spp*stressed*ts_length*search_effort 
                      +
                        f(comm_id, model = "iid",
                          hyper = list(prec = list(prior = "logtnormal",
                                                   param = c(0, 1)))) ,
                      data = inla_success_data2, 
                      family = "binomial",
                      lincomb = lincomb_summary,
                      control.compute = list(config = TRUE,cpo = TRUE),
                      control.predictor=list(compute=TRUE),
                      control.fixed = prior.fixed)

any(na.omit(inla_success2$cpo$failure)>=1)

fitted_rows_success2 <- which(is.na(inla_success_data2$sim_id))

marg_summary2 <- inla_success2$marginals.lincomb.derived[grepl("lc",names(inla_success2$marginals.lincomb.derived))] #subset marginals (equivalet to posteriors) to the linear combinations

summary_post2 <- lapply(marg_summary2,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x))
}) |>
  data.table::rbindlist() |>
  cbind(
    newdata[,c("n_spp","stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws)) #merge with newdata
  ) |>
  dplyr::mutate(ts_length = ts_length*100,
                stressed = ifelse(stressed == 1,"Yes","No")) #convert back to true ts length for visualisations

summary_range2 <- summary_post2 |>
  dplyr::group_by(n_spp,stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

#######################
# Fit multiAR
#######################

inla_success_data3 <- subset(summary_data,metric == "multiAR") %>%
  dplyr::bind_rows(pred_success_data)

inla_success3 <- inla(threshold_crossed ~ n_spp*stressed*ts_length*search_effort 
                      +
                        f(comm_id, model = "iid",
                          hyper = list(prec = list(prior = "logtnormal",
                                                   param = c(0, 1)))) ,
                      data = inla_success_data3, 
                      family = "binomial",
                      lincomb = lincomb_summary,
                      control.compute = list(config = TRUE,cpo = TRUE),
                      control.predictor=list(compute=TRUE),
                      control.fixed = prior.fixed)

any(na.omit(inla_success3$cpo$failure)>=1)

fitted_rows_success3 <- which(is.na(inla_success_data3$sim_id))

marg_summary3 <- inla_success3$marginals.lincomb.derived[grepl("lc",names(inla_success3$marginals.lincomb.derived))] #subset marginals (equivalet to posteriors) to the linear combinations

summary_post3 <- lapply(marg_summary3,function(x){
  data.frame(".draw" = 1:draws,
             ".value" = inla.rmarginal(draws, x))
}) |>
  data.table::rbindlist() |>
  cbind(
    newdata[,c("n_spp","stressed","ts_length","search_effort")] |>
      dplyr::slice(rep(1:dplyr::n(), each = draws)) #merge with newdata
  ) |>
  dplyr::mutate(ts_length = ts_length*100,
                stressed = ifelse(stressed == 1,"Yes","No")) #convert back to true ts length for visualisations

summary_range3 <- summary_post3 |>
  dplyr::group_by(n_spp,stressed,ts_length,search_effort) |>
  tidybayes::median_qi(.width = c(.95, .8, .5))

#######################
# Save posteriors
#######################
saveRDS(inla_success1,"Results/models/motif1_threshold_multiJI.rds")
saveRDS(inla_success2,"Results/models/motif1_threshold_uniJI.rds")
saveRDS(inla_success3,"Results/models/motif1_threshold_multiAR.rds")

save(summary_post1,summary_post2,summary_post3,
     file = "Results/models/motif1_threshold_posteriors.RData")

save(summary_range1,summary_range2,summary_range3,
     file = "Results/models/motif1_threshold_ranges.RData")
