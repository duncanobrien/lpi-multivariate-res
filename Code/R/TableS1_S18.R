########################################################################################################################
# Supplementary Tables #
########################################################################################################################
# models may have to be loaded and coefficients extracted individually depending on RAM
require(dplyr)
require(tibble)
source("Code/R/coef_inla_fn.R")

############ 
# 5 spp models (Tables S1-S5)
############ 
inla_5_1 <- readRDS("Results/models/motif1_5_invasive_multiJI_model.rds")
inla_5_1_coefs <- coefs_inla(inla_5_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_1)
write.csv(inla_5_1_coefs,
          file = "Results/figures/supplementary_tables/tableS1.csv", 
          row.names = FALSE)

inla_5_2 <- readRDS("Results/models/motif1_5_invasive_uniJI_model.rds")
inla_5_2_coefs <- coefs_inla(inla_5_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_2)
write.csv(inla_5_2_coefs,
          file = "Results/figures/supplementary_tables/tableS2.csv", 
          row.names = FALSE)

inla_5_3 <- readRDS("Results/models/motif1_5_invasive_FI_model.rds")
inla_5_3_coefs <- coefs_inla(inla_5_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_3)
write.csv(inla_5_3_coefs,
          file = "Results/figures/supplementary_tables/tableS3.csv", 
          row.names = FALSE)

inla_5_4 <- readRDS("Results/models/motif1_5_invasive_MVI_model.rds")
inla_5_4_coefs <- coefs_inla(inla_5_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_4)
write.csv(inla_5_4_coefs,
          file = "Results/figures/supplementary_tables/tableS4.csv", 
          row.names = FALSE)

inla_5_5 <- readRDS("Results/models/motif1_5_invasive_multiAR_model.rds")
inla_5_5_coefs <- coefs_inla(inla_5_5)|>
  tibble::add_column(metric = c("multiAR",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_5)
write.csv(inla_5_5_coefs,
          file = "Results/figures/supplementary_tables/tableS5.csv", 
          row.names = FALSE)

############ 
# 15 spp models (Tables S6-S10)
############ 
inla_15_1 <- readRDS("Results/models/motif1_15_invasive_multiJI_model.rds")
inla_15_1_coefs <- coefs_inla(inla_15_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_1)
write.csv(inla_15_1_coefs,
          file = "Results/figures/supplementary_tables/tableS6.csv", 
          row.names = FALSE)

inla_15_2 <- readRDS("Results/models/motif1_15_invasive_uniJI_model.rds")
inla_15_2_coefs <- coefs_inla(inla_15_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_2)
write.csv(inla_15_2_coefs,
          file = "Results/figures/supplementary_tables/tableS7.csv", 
          row.names = FALSE)

inla_15_3 <- readRDS("Results/models/motif1_15_invasive_FI_model.rds")
inla_15_3_coefs <- coefs_inla(inla_15_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_3)
write.csv(inla_15_3_coefs,
          file = "Results/figures/supplementary_tables/tableS8.csv", 
          row.names = FALSE)

inla_15_4 <- readRDS("Results/models/motif1_15_invasive_MVI_model.rds")
inla_15_4_coefs <- coefs_inla(inla_15_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_4)
write.csv(inla_15_4_coefs,
          file = "Results/figures/supplementary_tables/tableS9.csv", 
          row.names = FALSE)

inla_15_5 <- readRDS("Results/models/motif1_15_invasive_multiAR_model.rds")
inla_15_5_coefs <- coefs_inla(inla_15_5)|>
  tibble::add_column(metric = c("multiAR",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_5)
write.csv(inla_15_5_coefs,
          file = "Results/figures/supplementary_tables/tableS10.csv", 
          row.names = FALSE)

############ 
# 25 spp models (Tables S11-S15)
############ 
inla_25_1 <- readRDS("Results/models/motif1_25_invasive_multiJI_model.rds")
inla_25_1_coefs <- coefs_inla(inla_25_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_1)
write.csv(inla_25_1_coefs,
          file = "Results/figures/supplementary_tables/tableS11.csv", 
          row.names = FALSE)

inla_25_2 <- readRDS("Results/models/motif1_25_invasive_uniJI_model.rds")
inla_25_2_coefs <- coefs_inla(inla_25_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_2)
write.csv(inla_25_2_coefs,
          file = "Results/figures/supplementary_tables/tableS12.csv", 
          row.names = FALSE)

inla_25_3 <- readRDS("Results/models/motif1_25_invasive_FI_model.rds")
inla_25_3_coefs <- coefs_inla(inla_25_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_3)
write.csv(inla_25_3_coefs,
          file = "Results/figures/supplementary_tables/tableS13.csv", 
          row.names = FALSE)

inla_25_4 <- readRDS("Results/models/motif1_25_invasive_MVI_model.rds")
inla_25_4_coefs <- coefs_inla(inla_25_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_4)
write.csv(inla_25_4_coefs,
          file = "Results/figures/supplementary_tables/tableS14.csv", 
          row.names = FALSE)

inla_25_5 <- readRDS("Results/models/motif1_25_invasive_multiAR_model.rds")
inla_25_5_coefs <- coefs_inla(inla_25_5)|>
  tibble::add_column(metric = c("multiAR",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_5)
write.csv(inla_25_5_coefs,
          file = "Results/figures/supplementary_tables/tableS15.csv", 
          row.names = FALSE)

############ 
# Threshold models (Tables S16-S18)
############ 
inla_success1 <- readRDS("Results/models/motif1_threshold_multiJI.rds")
inla_success2 <- readRDS("Results/models/motif1_threshold_uniJI.rds")
inla_success3 <- readRDS("Results/models/motif1_threshold_multiAR.rds")

inla_success1_coefs <- coefs_inla(inla_success1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,18)),.before = 1)
rm(inla_success1)
write.csv(inla_success1_coefs,
          file = "Results/figures/supplementary_tables/tableS16.csv", 
          row.names = FALSE)

inla_success2_coefs <- coefs_inla(inla_success2) |>
  tibble::add_column(metric = c("max_uniJI",rep(NA,18)),.before = 1)
rm(inla_success2)
write.csv(inla_success2_coefs,
          file = "Results/figures/supplementary_tables/tableS17.csv", 
          row.names = FALSE)

inla_success3_coefs <- coefs_inla(inla_success3) |>
  tibble::add_column(metric = c("max_uniJI",rep(NA,18)),.before = 1)
rm(inla_success3)
write.csv(inla_success3_coefs,
          file = "Results/figures/supplementary_tables/tableS18.csv", 
          row.names = FALSE)
