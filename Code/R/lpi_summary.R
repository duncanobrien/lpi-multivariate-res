require(data.table)
require(magrittr)

lpi_raw_dt <- read.csv("Data/lpi_pops_20190927.csv") %>%
  as.data.table() %>%
  data.table::melt(., measure.vars = grep("^X",colnames(.)),
                   variable.name = "Year", value.name = "Value") %>%
  .[,Year := as.numeric(sub("X","",Year))] %>%
  data.table::setorder(.,ID) %>%
  data.table::setcolorder(.,c("ID","Binomial","System","Year","Value")) %>%
  .[,Value := as.numeric(Value)] %>%
  #.[,Value := ifelse(Value == "NULL",NA,Value)] %>%
  .[,Migratory := ifelse(is.na(Migratory),"NA",Migratory)] %>% 
  #sporadic NAs in migratory column. Convert to character to allow missing dates to 
  #be dropped efficiently (is the only additional column with NAs)
  na.omit(.) 

tot_summary_lpi <- copy(lpi_raw_dt) %>%
  .[,.(ts_len = length(Year),
       mean_val = mean(Value),
       year_range = paste(min(Year),max(Year),sep = "-")),
    by=c("ID","Binomial","System")] %>% #summarise total ts length, mean value and year range per ID
  .[,`:=`(
    system_mean_len = mean(ts_len),
    system_max_len = max(ts_len),
    system_min_len = min(ts_len)
  ),by="System"] %>% #summarise total ts length per system
  .[,`:=`(
    tot_mean_len = mean(ts_len),
    tot_max_len = max(ts_len),
    to_min_len = min(ts_len))]

hist(tot_summary_lpi$ts_len)
quantile(tot_summary_lpi$ts_len,c(0.05,0.25,0.5,0.75,0.95))

tt <- round(rnorm(100,mean=50,15))
kk <- sapply(tt,FUN = function(x){
  round(x*pnbinom(q=x*0.75,size = x, prob = 0.5))
})


  