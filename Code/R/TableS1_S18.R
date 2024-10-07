########################################################################################################################
# Supplementary Tables #
########################################################################################################################
# models may have to be loaded and coefficients extracted individually depending on RAM
require(dplyr)
require(tibble)
require(tinytable)
source("Code/R/auxiliary_functions/coef_inla_fn.R")

############ 
# 5 spp models (Tables S1-S5)
############ 
inla_5_1 <- readRDS("Results/models/motif1_5_invasive_multiJI_model.rds")
inla_5_1_coefs <- coefs_inla(inla_5_1) |>
  tibble::add_column(metric = c("multiJI",rep("NA",21)),.before = 1) |>
  tibble::add_column(n_spp = c("5",rep("NA",21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS1.docx", overwrite = TRUE)
rm(inla_5_1)

inla_5_2 <- readRDS("Results/models/motif1_5_invasive_uniJI_model.rds")
inla_5_2_coefs <- coefs_inla(inla_5_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1)|>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS2.docx", overwrite = TRUE)
rm(inla_5_2)

inla_5_3 <- readRDS("Results/models/motif1_5_invasive_FI_model.rds")
inla_5_3_coefs <- coefs_inla(inla_5_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS3.docx", overwrite = TRUE)
rm(inla_5_3)

inla_5_4 <- readRDS("Results/models/motif1_5_invasive_MVI_model.rds")
inla_5_4_coefs <- coefs_inla(inla_5_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS4.docx", overwrite = TRUE)
rm(inla_5_4)

inla_5_5 <- readRDS("Results/models/motif1_5_invasive_multiAR_model.rds")
inla_5_5_coefs <- coefs_inla(inla_5_5)|>
  tibble::add_column(metric = c("multiAR",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(5,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS5.docx", overwrite = TRUE)
rm(inla_5_5)

############ 
# 15 spp models (Tables S6-S10)
############ 
inla_15_1 <- readRDS("Results/models/motif1_15_invasive_multiJI_model.rds")
inla_15_1_coefs <- coefs_inla(inla_15_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS6.docx", overwrite = TRUE)
rm(inla_15_1)

inla_15_2 <- readRDS("Results/models/motif1_15_invasive_uniJI_model.rds")
inla_15_2_coefs <- coefs_inla(inla_15_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS7.docx", overwrite = TRUE)
rm(inla_15_2)

inla_15_3 <- readRDS("Results/models/motif1_15_invasive_FI_model.rds")
inla_15_3_coefs <- coefs_inla(inla_15_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS8.docx", overwrite = TRUE)
rm(inla_15_3) 

inla_15_4 <- readRDS("Results/models/motif1_15_invasive_MVI_model.rds")
inla_15_4_coefs <- coefs_inla(inla_15_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS9.docx", overwrite = TRUE)
rm(inla_15_4)

inla_15_5 <- readRDS("Results/models/motif1_15_invasive_multiAR_model.rds")
inla_15_5_coefs <- coefs_inla(inla_15_5)|>
  tibble::add_column(metric = c("multiAR",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(15,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS10.docx", overwrite = TRUE)
rm(inla_15_5)

############ 
# 25 spp models (Tables S11-S15)
############ 
inla_25_1 <- readRDS("Results/models/motif1_25_invasive_multiJI_model.rds")
inla_25_1_coefs <- coefs_inla(inla_25_1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS11.docx", overwrite = TRUE)
rm(inla_25_1)

inla_25_2 <- readRDS("Results/models/motif1_25_invasive_uniJI_model.rds")
inla_25_2_coefs <- coefs_inla(inla_25_2)|>
  tibble::add_column(metric = c("mean_uniJI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS12.docx", overwrite = TRUE)
rm(inla_25_2)

inla_25_3 <- readRDS("Results/models/motif1_25_invasive_FI_model.rds")
inla_25_3_coefs <- coefs_inla(inla_25_3)|>
  tibble::add_column(metric = c("FI",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS13.docx", overwrite = TRUE)
rm(inla_25_3)

inla_25_4 <- readRDS("Results/models/motif1_25_invasive_MVI_model.rds")
inla_25_4_coefs <- coefs_inla(inla_25_4)|>
  tibble::add_column(metric = c("mvi",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS14.docx", overwrite = TRUE)
rm(inla_25_4)

inla_25_5 <- readRDS("Results/models/motif1_25_invasive_multiAR_model.rds")
inla_25_5_coefs <- coefs_inla(inla_25_5)|>
  tibble::add_column(metric = c("multiAR",rep(NA,21)),.before = 1) |>
  tibble::add_column(n_spp = c(25,rep(NA,21)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = c(1,2), rowspan = 22,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 3, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 2, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS15.docx", overwrite = TRUE)
rm(inla_25_5)

############ 
# Threshold models (Tables S16-S18)
############ 
inla_success1 <- readRDS("Results/models/motif1_threshold_multiJI.rds")
inla_success2 <- readRDS("Results/models/motif1_threshold_uniJI.rds")
inla_success3 <- readRDS("Results/models/motif1_threshold_multiAR.rds")

inla_success1_coefs <- coefs_inla(inla_success1) |>
  tibble::add_column(metric = c("multiJI",rep(NA,18)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = 1, rowspan = 19,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 2, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 1, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS16.docx", overwrite = TRUE)
rm(inla_success1)

inla_success2_coefs <- coefs_inla(inla_success2) |>
  tibble::add_column(metric = c("max_uniJI",rep(NA,18)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = 1, rowspan = 19,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 2, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 1, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS17.docx", overwrite = TRUE)
rm(inla_success2)

inla_success3_coefs <- coefs_inla(inla_success3) |>
  tibble::add_column(metric = c("max_uniJI",rep(NA,18)),.before = 1) |>
  tinytable::tt(theme = "grid") |>
  tinytable::style_tt(i = 1, j = 1, rowspan = 19,alignv = "m") |>
  tinytable::style_tt(i = c(1), j = 2, colspan = 5,align = "l") |>
  tinytable::style_tt(i = c(18), j = 1, colspan = 6,align = "l") |>
  tinytable::save_tt("Results/figures/supplementary_tables/tableS18.docx", overwrite = TRUE)
rm(inla_success3)
