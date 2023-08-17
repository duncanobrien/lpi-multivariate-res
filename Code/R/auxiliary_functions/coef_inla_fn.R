#' Extract coefficients of an INLA model

#' @param obj An INLA model object
#' @returns A data.frame of model coefficients
#' 
coefs_inla <- function(obj){
  summ <- summary(obj)
  
  fixed <- summ$fixed |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Parameter") |>
    dplyr::select("Parameter","0.5quant","0.025quant","0.975quant","kld") |>
    `colnames<-`(c("Parameter","Median","CI_low","CI_high","kld")) 
  
  random <- summ$hyperpar
  random <- random[!grepl("Gaussian observation",rownames(random)),] |> #drop unnecessary estimate
    tibble::rownames_to_column(var = "Parameter") |>
    dplyr::select("Parameter","0.5quant","0.025quant","0.975quant") |>
    `colnames<-`(c("Parameter","Median","CI_low","CI_high")) |>
    dplyr::mutate(kld = NA)
  
  out <- rbind(rbind(c("Fixed effects",rep(NA,4)),fixed),
               rbind(c("Random effects",rep(NA,4)),random))
  
  return(out)
}
