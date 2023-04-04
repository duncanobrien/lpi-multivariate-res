require(tidyverse)
require(data.table)
require(brms)

test_summary_data <- data.frame(sim_id = unite(expand.grid(seq(1:5),seq(1:5),seq(1:5)),sim_id)) |>
  slice(rep(1:n(),each = 4)) %>%
  mutate(metric = ,
         trend =
           ,corrected_trend =
           ,threshold_crossed =
           ,target_node =
           ,jac_collapse =
           ,model,stressed =
           ,ts_length = 
           ,search_effort = rep(c(0.25,0.5,0.75,0.95),each = 125))

test_summary_data <- expand.grid(sim_id = unite(expand.grid(seq(1:5),seq(1:5),seq(1:5)),sim_id)[[1]],
            search_effort = c(0.25,0.5,0.75,0.95),
            model = c("harvest","invasive"),
            stressed = c(1,0),
            ts_length = c(3,9,21,45,70),
            metric = c("multiJI","mean_uniJI","max_uniJI","FI","mvi")) %>%
  rowwise() %>%
  mutate(corrected_trend = ifelse(metric %in% c("multiJI","max_uniJI"), rnorm(1,0.5,0.25), rnorm(1,0.1,0.4)),
         threshold_crossed = ifelse(metric %in% c("multiJI","max_uniJI"), rbinom(1,1,0.95),
                                    ifelse(metric %in% c("mean_uniJI"), rbinom(1,1,0.45),NA))) %>%
  as.data.table() %>%
  .[,success := ifelse(threshold_crossed == 1 & stressed ==1, 1,
                       ifelse(threshold_crossed == 0 & stressed == 0, 1,
                              ifelse(is.na(threshold_crossed), NA,0)))] %>%
  .[, c("motif", "community","sim") := data.table::tstrsplit(sim_id, "_", fixed=TRUE)] %>%
  as.data.frame()

its <- 2000
thn <- 0.0005*its
wrmup <- 0.1*its

################################################
# Model threshold
################################################
get_prior(brms::bf(success ~ metric*ts_length*search_effort*model - 1 + (1|motif/community)),
          data = test_summary_data,
          family = binomial(link = "logit"))

mod_summary_data <- as.data.table(test_summary_data) %>%
  .[,.(total_success = sum(as.numeric(success)),
        trials = length(unique(sim_id)),by = c("metric","ts_length","search_effort","model"))] %>%
  as.data.frame()

success_mod <- brms::brm(brms::bf(success ~ metric*ts_length*search_effort*model - 1 + (1|motif/community)), 
                         data = test_summary_data,
                         iter = its,
                         thin = thn,
                         warmup = wrmup,
                         prior= c(prior(normal(0, 1.2), class = b),
                                  prior(normal(1, 1),class = sd)),
                         family = bernoulli(link = "logit"), 
                         chains = 2,
                         control = list(adapt_delta = .975, max_treedepth = 20),
                         seed = 12345, cores = 2,sample_prior = TRUE)

trend_mod <- brms::brm(brms::bf(corrected_trend ~ metric*ts_len*search_effort*model1 + (1|motif/community)), 
                       data = test_summary_data,
                       iter = its,
                       thin = thn,
                       warmup = wrmup,
                       prior= c(prior(normal(0, 2), class = b),
                                prior(exp(1),class = sd)),
                       family = gaussian(), 
                       chains = 4,
                       control = list(adapt_delta = .975, max_treedepth = 20),
                       seed = 12345, cores = 4,sample_prior = TRUE)
