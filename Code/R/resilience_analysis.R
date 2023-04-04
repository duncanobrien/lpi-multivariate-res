########################################################################
## Preamble ##
########################################################################
require(brms)
require(data.table)
require(magrittr)

################################################
# Load resilience data
################################################
full_path <- paste("Data/resilience/full") #define stressed model path
full_files <- fs::dir_ls(path = full_path, recurse = TRUE,glob = "*.csv")

summary_path <- paste("Data/resilience/summary") #define stressed model path
summary_files <- fs::dir_ls(path = summary_path, recurse = TRUE,glob = "*.csv")

full_data <- vroom::vroom(full_files, id = "sample")
summary_data <- vroom::vroom(summary_files, id = "sample")

full_data <- full_data %>%
  melt(.,measure.vars = c("multiJI","mean_uniJI","max_uniJI","FI","mvi"),
       variable.name = "metric",value.name = "metric_value") %>% #pivot resilience metrics longer
  .[, c("motif", "community","sim") := data.table::tstrsplit(sim_id, "_", fixed=TRUE)]

summary_data <- summary_data %>%
  .[,success := ifelse(threshold_crossed == 1 & stressed ==1, 1,
                       ifelse(threshold_crossed == 0 & stressed == 0, 1,
                              is.na(threshold_crossed), NA,0))] %>%
  .[, c("motif", "community","sim") := data.table::tstrsplit(sim_id, "_", fixed=TRUE)]

its <- 10000
thn <- 0.0005*its
wrmup <- 0.1*its

################################################
# Model threshold
################################################

success_mod <- brms::brm(brms::bf(success ~ metric*ts_len*search_effort*model - 1 + (1|motif/community)), 
                                    data = summary_data,
                                    iter = its,
                                    thin = thn,
                                    warmup = wrmup,
                                    prior= c(prior(normal(0, 1.2), class = b),
                                             prior(normal(1, 1),class = sd)),
                                    family = binomial(link = "logit"), 
                                    chains = 4,
                                    control = list(adapt_delta = .975, max_treedepth = 20),
                                    seed = 12345, cores = 4,sample_prior = TRUE)

################################################
# Model trend
################################################
trend_mod <- brms::brm(brms::bf(metric_value ~ metric*ts_len*search_effort*model*time - 1 + (1|motif/community)), 
                         data = subset(full_data,metric %in% "multiJI","mean_uniJI","FI","mvi"),
                         iter = its,
                         thin = thn,
                         warmup = wrmup,
                         prior= c(prior(normal(0, 2), class = b),
                                  prior(exp(1),class = sd)),
                         family = gaussian(), 
                         chains = 4,
                         control = list(adapt_delta = .975, max_treedepth = 20),
                         seed = 12345, cores = 4,sample_prior = TRUE)

