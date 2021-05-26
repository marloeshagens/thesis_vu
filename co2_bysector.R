setwd("/cloud/project")
cat("\014")        
rm(list = ls())

####### install packages
#######
# install.packages("Synth")
# install.packages("ArCo")
# install.packages("xlsx")
# install.packages("ggplot2")
# install.packages("moments")
# install.packages("dplyr")
# install.packages("tseries") # time series
# install.packages("urca") # contains adf unit root (ur.df) & kpss trend stationarity tests (ur.kpss)
# install.packages("reshape2")
# install.packages("Java")
# install.packages("tidyverse") #includes readxl function
# install.packages(glmnet) #contains LASSO function (Friedman et al., 2010).

library(xlsx)
library(ggplot2)
library(moments)
library(tseries) # time series
library(urca) # contains adf unit root (ur.df) & kpss trend stationarity tests (ur.kpss)
library(reshape2) # contains acast function, to reshape df
library(Synth) # Synthetic Control method
library(ArCo) # Artificial Counterfactual method
library(readxl) # from tidyverse
# library(tidyverse) #complete tidyverse, not necessary
library(dplyr)
library(glmnet) #contains LASSO function (Friedman et al., 2010).

alpha = 0.05

#######
###### ArCo on UK CO2 emissions 1990-2018; t0=2013

# data source: CAIT (2020) - extracted 19 may 2020 
# sector: Electricity/Heat
# gas/emissions: CO2
# measurement: MtCO_2_e

df <- read.csv("historical_emissions.csv", header = TRUE)

# data format for ArCo
rownames(df) <- df[,1]
df <- df[-c(1:5)]
for (col in 1:ncol(df)){
  colnames(df)[col] <-  sub("X","", colnames(df)[col])
}
tdf <- t(df) ##-- i.e.,  a[i, j] == ta[j, i] for all i,j 
for(j in seq(ncol(df)))
  if(! all(df[, j] == tdf[j, ])) stop("wrong transpose")
tdf<- as.data.frame(tdf)
tdf<- tdf[c(29:1),]
#tdf<- as.ts(tdf, start = "1990", end = "2018")

# UK with all 27 EU member states + Norway (EEA-EFTA state) 
# excl. Iceland (EEA-EFTA state) & Liechtenstein (because of NA and 0 values)
data <- subset(tdf,select=
                c( "United Kingdom","Austria","Belgium","Bulgaria","Croatia",
                   "Cyprus","Czech Republic","Denmark","Estonia","Finland",
                   "France","Germany","Greece","Hungary","Ireland", "Italy",
                   "Latvia","Lithuania","Luxembourg","Malta","Netherlands", 
                   "Poland","Portugal","Romania","Slovakia","Slovenia",
                   "Spain", "Sweden")) # no "Norway" yet, maybe add? 
# data <- subset(data[,-4]) #remove bulgaria
input <- list(data)

t0 = 24 # 2013 is on row 24

# ArCo first stage is fitted using LASSO (cv.glmnet, predict)
# glmnet automatically normalizes data & transforms data back (normalization is needed for lasso)
arco <- fitArCo(data = input,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = 1, t0 = t0, boot.cf=TRUE,R=200,l=3 
                #VCOV.type = "nw", prewhitening.kernel = TRUE
                )
plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "UK ArCo", ylab = "CO2 emissions (MtCO₂e)")
delta<-arco$delta 
p<-arco$p.value 

# UNIT ROOT & TREND STATIONARITY TESTS
adf <- lapply(data, adf.test)
kpss <- lapply(data, kpss.test, null = "Trend")
unitroot_table<- function(df){
  p <- ncol(df)
  df_stats <- data.frame(var=names(df),
                         adf.pvalue=sapply(df, function(v) adf.test(ts(v),alternative = "stationary")$p.value),
                         kpss.pvalue=sapply(df, function(v) kpss.test(ts(v),null = "Trend")$p.value)
                         )
  df_stats$unitroot_adf <- ifelse(df_stats$adf.pvalue < 0.05, 0, 1)
  df_stats$unitroot_kpss <- ifelse(df_stats$kpss.pvalue < 0.05, 1, 0)
  df_stats$conclusion <- ifelse(df_stats$unitroot_kpss + df_stats$unitroot_adf > 1, 
                              "contains unit root", "trend stationary or unclear")
  row.names(df_stats) <- c()
  df_stats
}
tests <- unitroot_table(data) # luxembourg and sweden are trend-stationary

#difference the data
diff_data <- diff(as.matrix(data))
diff_tests <- unitroot_table(as.data.frame(diff_data)) # only bulgaria still contains unit root
diff_data <- subset(diff_data[,-4]) #remove bulgaria

#apply ArCo on  differenced data
input2 <- list(diff_data)
arco2 <- fitArCo(data = input2,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = 1, t0 = 24, boot.cf=TRUE, R=200, l=3, 
                #VCOV.type = "nw", prewhitening.kernel = TRUE
                )
plot(arco2, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "UK ArCo", ylab = "diff. CO2 emissions (MtCO₂e)")
delta2<-arco2$delta 
p2<-arco2$p.value 
p2<0.05 #TRUE 

# mimicking Leroutier's simulation
# check the reason for the exclusion of some countries
data_leroutier <- subset(tdf, 
                select=c("United Kingdom","Austria","Belgium",
                         "Czech Republic","Denmark","Finland","France",
                         "Hungary","Ireland", "Italy","Netherlands",
                         "Poland","Portugal","Slovakia","Spain","Sweden")) 
diff_data_leroutier <- subset(diff_data, 
                              select=c("United Kingdom","Austria","Belgium",
                                       "Czech Republic","Denmark","Finland","France",
                                       "Hungary","Ireland", "Italy","Netherlands",
                                       "Poland","Portugal","Slovakia","Spain","Sweden")) 
input3 <- list(data_leroutier) 
arco3 <- fitArCo(data = input3, fn = cv.glmnet, p.fn = predict, treated.unit = 1 , t0 = t0)
plot(arco3, display.fitted=TRUE, alpha = 0.05, main= "UK ArCo - same countries as Leroutier", ylab = "CO2 emissions (MtCO₂e)")
delta3<-arco3$delta 
p3<-arco3$p.value 
p3<alpha #TRUE 

input4 <- list(diff_data_leroutier) 
arco4 <- fitArCo(data = input4, fn = cv.glmnet, p.fn = predict, treated.unit = 1 , t0 = t0)
plot(arco4, display.fitted=TRUE, alpha = 0.05, main= "UK ArCo - same countries as Leroutier", ylab = "diff. CO2 emissions (MtCO₂e)")
delta4<-arco4$delta 
p4<-arco4$p.value 
p4<alpha #TRUE 

### PLOTS 
# # get data out of arco (see plot.arcofit.R file)
# Y = matrix(arco$data[[1]][, 1], ncol = 1)
# Y_fit = arco$fitted.values
# cf=arco$cf
# boot.cf=arco$boot.cf
# confidence.bands=TRUE
# display.fitted=TRUE
# alpha = 0.05
# y.min=apply(rbind(Y,cf),2,min,na.rm=TRUE) 
# y.max=apply(rbind(Y,cf),2,max,na.rm=TRUE)
# plot(Y, type = "l", xlab = "Time",ylim=c(y.min,y.max))
#   


#######
# data information: 
# ELECTRICITY & HEAT (CAIT, 2020)
# co2 emissions per sector/country/year 
# Variable time span: 1990-2018
# extraction date: 19-05-2021
# extraction link: https://www.climatewatchdata.org/data-explorer/historical-emissions?historical-emissions-data-sources=cait&historical-emissions-end_year=2018&historical-emissions-gases=All%20Selected%2Cco2&historical-emissions-regions=GBR%2CEUU&historical-emissions-sectors=electricity-heat&historical-emissions-start_year=1990&page=1&sort_col=country&sort_dir=ASC#data
# Data published by	CAIT Climate Data Explorer via. Climate Watch
# Carbon dioxide (CO₂) emissions only given for Electricity/Heat sector
# measured in MtCO₂e per year.
# This data is published by country and sector from the CAIT Climate Data Explorer, and downloaded from the Climate Watch Portal. Available here: https://www.climatewatchdata.org/data-explorer/historical-emissions

# General Link to the Data Explorer:	https://www.climatewatchdata.org/data-explorer/historical-emissions



# ####### other data option; same source, but only 1990-2016. 
# # deselect below
# 
# df <- read.csv("co-emissions-by-sector.csv", header = TRUE)
# 
# # select only electricity and heat sector
# df2<-dcast(df, Year ~ Entity, value.var="Electricity...Heat..CAIT..2020.")
# 
# # UK with all 27 EU member states + Norway (EEA-EFTA state) 
# # excl. Iceland (EEA-EFTA state) & Liechtenstein (because of NA values and 0 values)
# df3 <- subset(df2,select=
#                 c( "United Kingdom","Austria","Belgium","Bulgaria","Croatia",
#                    "Cyprus","Czechia","Denmark","Estonia","Finland",
#                    "France","Germany","Greece","Hungary","Ireland", "Italy",
#                    "Latvia","Lithuania","Luxembourg","Malta","Netherlands", 
#                    "Norway","Poland","Portugal","Romania","Slovakia","Slovenia",
#                    "Spain", "Sweden"))
# data <- list(df3)
# 
# t0 <- which(grepl(2013, df2$Year))
# 
# # ArCo first stage is fitted using LASSO (cv.glmnet, predict)
# arco <- fitArCo(data = data,  fn = cv.glmnet, p.fn = predict, treated.unit = 1 , t0 = t0, 
#                 boot.cf=TRUE,R=200,l=3, VCOV.type = "nw", prewhitening.kernel = TRUE)
# plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05)
# delta<-arco$delta 
# p<-arco$p.value 
# 
# # ELECTRICITY & HEAT (CAIT, 2020)
# # extracted from https://ourworldindata.org/emissions-by-sector#sector-by-sector-where-do-global-greenhouse-gas-emissions-come-from
# # Variable time span	1990 – 2016
# # Unit conversion factor for chart	1000000
# # Data published by	CAIT Climate Data Explorer via. Climate Watch
# # Link	https://www.climatewatchdata.org/data-explorer/historical-emissions
# # Carbon dioxide (CO₂) emissions broken down by sector, measured in tonnes per year. Further information on sector definitions is available here: https://ourworldindata.org/ghg-emissions-by-sector
# # This data is published by country and sector from the CAIT Climate Data Explorer, and downloaded from the Climate Watch Portal. Available here: https://www.climatewatchdata.org/data-explorer/historical-emissions
# # check for more info: https://ourworldindata.org/ghg-emissions-by-sector 
