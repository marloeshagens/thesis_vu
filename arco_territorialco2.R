### DIFFERENT DATA - ARCO WORKS

# setwd("/Users/marloes/Documents/thesis Climate Econometrics R project")
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

#######
#Begin
CPS_year <- 2013
UK <- "United.Kingdom"

#df <- read_xlsx("MH_National_Carbon_Emissions_2020v1.0.xlsx", 3)
emissions_data <- read.xlsx("MH_National_Carbon_Emissions_2020v1.0.xlsx", 3)
t0 <- which(grepl(CPS_year, emissions_data$Year)) # row: year of implementation
treated <- which(colnames(emissions_data)==UK)

# UK with all 27 EU member states + Iceland and Norway (EEA-EFTA states)
# excl Liechtenstein (because of NA values)
df <- subset(emissions_data, 
             select=c("United.Kingdom","Austria","Belgium","Bulgaria","Croatia",
                      "Cyprus","Czech.Republic","Denmark","Estonia","Finland",
                      "France","Germany","Greece","Hungary","Ireland", "Italy",
                      "Latvia","Lithuania","Luxembourg","Malta","Netherlands", 
                      "Poland","Portugal","Romania","Slovakia","Slovenia",
                      "Spain", "Sweden",
                      "Iceland","Norway"))
data <- list(df)

# ArCo can be fitted using linear regression functions below:
fn=function(X,y){
  return(lm(y~X))
}
p.fn=function(model,newdata){
  b=coef(model)
  return(cbind(1,newdata) %*% b)}

arco1 <- fitArCo(data = data, fn = fn, p.fn = p.fn, treated.unit = 1 , t0 = t0)
plot(arco1)
plot(arco1, display.fitted=TRUE, alpha = 0.05)
delta<-arco1$delta 
p<-arco1$p.value 

#using LASSO 
arco2 <- fitArCo(data = data,  fn = cv.glmnet, p.fn = predict, treated.unit = 1 , t0 = t0, 
                 boot.cf=TRUE,R=200,l=3)
                  #VCOV.type = "nw", prewhitening.kernel = TRUE,
plot(arco2, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05)

# same selected countries as Leroutier 
# check the reason for the exclusion of some countries
df_ML <- subset(emissions_data, 
             select=c("United.Kingdom","Austria","Belgium",
                      "Czech.Republic","Denmark","Finland","France",
                      "Hungary","Ireland", "Italy","Netherlands",
                      "Poland","Portugal","Slovakia","Sweden","Switzerland")) 
data_ML <- list(df_ML) 
arco3 <- fitArCo(data = data_ML, fn = fn, p.fn = p.fn, treated.unit = 1 , t0 = t0)
plot(arco3, display.fitted=TRUE, alpha = 0.05)
delta3<-arco3$delta 
p3<-arco3$p.value 

arco4 <- fitArCo(data = data_ML,  fn = cv.glmnet, p.fn = predict, treated.unit = 1, 
                 t0 = t0, 
                 boot.cf=TRUE,R=200,l=1)
                  #VCOV.type = "nw", prewhitening.kernel = TRUE,
plot(arco4, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05)

#### INFORMATION ON DATA 
#########################################################################
# Fossil CO2 emissions by country (territorial)
# All values in million tonnes of carbon per year. For values in million tonnes of CO2 per year, multiply the values below by 3.664
# 1MtC = 1 million tonne of carbon = 3.664 million tonnes of CO2
# Cite as: Friedlingstein et al. 2020
# Methods: Full details of the method are described in Friedlingstein et al (2020).
# (1) National estimates include emissions from fossil fuel combustion and oxidation and cement production and excludes emissions from bunker fuels. World totals include emissions from bunker fuels. 
# (2) Bunker fuels: Emissions from fuels used for international aviation and maritime transport
# (3) The disaggregations of regions (e.g. the former Soviet Union prior to 1992) are based on the shares of emissions in the first year after the countries are disaggregated (e.g., 1992 for the Former Soviet Union).
# (4) The statistical difference presented on column HX is the difference between the world emissions and the sum of the emissions for each countries and for the bunker fuels.