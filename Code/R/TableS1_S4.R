########################################################################################################################
# Supplementary Tables #
########################################################################################################################
# models may have to be loaded and coefficients extracted individually depending on RAM
require(dplyr)
require(tibble)
source("Code/R/coef_inla_fn.R")

############ 
# Table S1
############ 
inla_5_1 <- readRDS("Results/models/motif1_5_invasive_multiJI_model.rds")
inla_5_1_coefs <- coefs_inla(inla_5_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_1)

inla_5_2 <- readRDS("Results/models/motif1_5_invasive_uniJI_model.rds")
inla_5_2_coefs <- coefs_inla(inla_5_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_2)

inla_5_3 <- readRDS("Results/models/motif1_5_invasive_FI_model.rds")
inla_5_3_coefs <- coefs_inla(inla_5_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_3)

inla_5_4 <- readRDS("Results/models/motif1_5_invasive_MVI_model.rds")
inla_5_4_coefs <- coefs_inla(inla_5_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)
rm(inla_5_4)

write.csv(rbind(inla_5_1_coefs,inla_5_2_coefs,inla_5_3_coefs,inla_5_4_coefs),
          file = "Results/figures/supplementary_tables/tableS1.csv", 
          row.names = FALSE)

############ 
# Table S2
############ 
inla_15_1 <- readRDS("Results/models/motif1_15_invasive_multiJI_model.rds")
inla_15_1_coefs <- coefs_inla(inla_15_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_1)

inla_15_2 <- readRDS("Results/models/motif1_15_invasive_uniJI_model.rds")
inla_15_2_coefs <- coefs_inla(inla_15_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_2)

inla_15_3 <- readRDS("Results/models/motif1_15_invasive_FI_model.rds")
inla_15_3_coefs <- coefs_inla(inla_15_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_3)

inla_15_4 <- readRDS("Results/models/motif1_15_invasive_MVI_model.rds")
inla_15_4_coefs <- coefs_inla(inla_15_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1)
rm(inla_15_4)

write.csv(rbind(inla_15_1_coefs,inla_15_2_coefs,inla_15_3_coefs,inla_15_4_coefs),
          file = "Results/figures/supplementary_tables/tableS2.csv", 
          row.names = FALSE)

############ 
# Table S3
############ 
inla_25_1 <- readRDS("Results/models/motif1_25_invasive_multiJI_model.rds")
inla_25_1_coefs <- coefs_inla(inla_25_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_1)

inla_25_2 <- readRDS("Results/models/motif1_25_invasive_uniJI_model.rds")
inla_25_2_coefs <- coefs_inla(inla_25_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_2)

inla_25_3 <- readRDS("Results/models/motif1_25_invasive_FI_model.rds")
inla_25_3_coefs <- coefs_inla(inla_25_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_3)

inla_25_4 <- readRDS("Results/models/motif1_25_invasive_MVI_model.rds")
inla_25_4_coefs <- coefs_inla(inla_25_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1)
rm(inla_25_4)

write.csv(rbind(inla_25_1_coefs,inla_25_2_coefs,inla_25_3_coefs,inla_25_4_coefs),
          file = "Results/figures/supplementary_tables/tableS3.csv", 
          row.names = FALSE)

############ 
# Table S4
############ 
inla_success1 <- readRDS("Results/models/motif1_threshold_multiJI.rds")
inla_success2 <- readRDS("Results/models/motif1_threshold_uniJI.rds")

inla_success1_coefs <- coefs_inla(inla_success1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,18)),.before = 1)
rm(inla_success1)

inla_success2_coefs <- coefs_inla(inla_success2) |>
  tibble::add_column(metric = c("max_uniJI",rep(NA,18)),.before = 1)
rm(inla_success2)

write.csv(rbind(inla_success1_coefs,inla_success2_coefs),
          file = "Results/figures/supplementary_tables/tableS4.csv", 
          row.names = FALSE)
