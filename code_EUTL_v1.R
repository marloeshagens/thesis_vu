# setwd("/Users/marloes/Documents/thesis Climate Econometrics R project")
setwd("/cloud/project")
cat("\014")        
rm(list = ls())

################################################################################  
############################### 0. LOAD PACKAGES ############################### 
################################################################################ 
####### install packages #######
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

####### load packages #######
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

################################################################################  
############################ 1. LOAD & PREPARE DATA ############################
################################################################################  

####1.1. emissions data (EU)####

#select only the power installations
compliance_data <- read.csv("compliance.csv", header = TRUE)
installation_data <- read.csv("installation.csv", header = TRUE)
joined_df <- merge(compliance_data, installation_data, by.x = "installation_id", 
                   by.y = "id", all.x = TRUE, all.y = FALSE)
df <- subset(joined_df, activity_id == 1 | activity_id == 20) # only activity 1 and 20; power installations

#split & select United Kingdom and rest of EU data set 
df_GB <- subset(df, country_id == "GB")
df_exclGB <- subset(df, country_id != "GB")

#### 1.2. data for predictors used for the SC method ####

# 1.2.1 coal-to-gas price ratio

# coal prices

# gas prices

# 1.2.2 lignite resources 
# 1.2.3 Residual load per capita 
# 1.2.4 Emission from opted-out plants
# 1.2.5 Number of degree days
# 1.2.6 Per capita capacity for combustible fuels, gas and coal
# 1.2.7 Growth in per capita wind and solar capacity

################################################################################  
############################### 2. DESCRIPTIVES ################################ 
################################################################################  

#### descriptives and plots #######
## the whole EU
df_aggregate <- aggregate(df$verified, by=list(year=df$year), FUN=sum, na.rm=TRUE, na.action=NULL)
df_aggregate <- df_aggregate[-c(16:26), ]
# augmented dickey-fueller test - stationarity

plot(ts(df_aggregate$x, start = 2005 , end = 2019),
        ylab = "CO2 emissions (t)", xlab = "Year",
        type = "l", 
        main = "Aggregated Verified Emissions in the EU-ETS",
        ylim=c(800000000, 1600000000),
        xlim=c(2005,2019))

## Great Britain
df_GB_aggregate <- aggregate(df_GB$verified, by=list(year=df_GB$year), FUN=sum, na.rm=TRUE, na.action=NULL)
df_GB_aggregate <- df_GB_aggregate[-c(16), ]

plot(ts(df_GB_aggregate$x, start = 2005 , end = 2019),
     ylab = "CO2 emissions (t)", xlab = "Year",
     type = "l", 
     main = "Aggregated Verified Emissions in the GB, EU-ETS",
     #ylim=c(800000000, 1600000000),
     xlim=c(2005,2019))
abline(v=2013, col="grey", lty=2, lwd=2)

## rest of the EU
df_exclGB_aggregate <- aggregate(df_exclGB$verified, by=list(year=df_exclGB$year), FUN=sum, na.rm=TRUE, na.action=NULL)
df_exclGB_aggregate <- df_exclGB_aggregate[-c(16:26), ]

plot(ts(df_exclGB_aggregate$x, start = 2005 , end = 2019),
     ylab = "CO2 emissions (t)", xlab = "Year",
     type = "l", 
     main = "Aggregated Verified Emissions in the EU-ETS, excl GB",
     #ylim=c(800000000, 1600000000),
     xlim=c(2005,2019))

#### Normalize function to compare the data ####
# normalizing function, such we can compare the data to each
normalize <- function(x) {
        return ((x - min(x)) / (max(x) - min(x))) 
        }
df_GB_aggregate$x_norm <- normalize(df_GB_aggregate$x)
df_aggregate$x_norm <- normalize(df_aggregate$x)
df_exclGB_aggregate$x_norm <- normalize(df_exclGB_aggregate$x)

# normalized data plots in 1 
plot(ts(df_GB_aggregate$x_norm, start = 2005 , end = 2019),
     ylab = "Norm. CO2 emissions", 
     xlab = "Year",
     type = "l", 
     col = "red",
     main = "Normalized Verified Emissions Data in the EU-ETS",
     xlim=c(2005,2019))
lines(ts(df_aggregate$x_norm, start = 2005 , end = 2019), type = "l", col = "grey")
lines(ts(df_exclGB_aggregate$x_norm, start = 2005 , end = 2019), type = "l", col = "blue")
abline(v=2013, col="grey", lty=2, lwd=2)
legend("topright",
       legend=c(expression(EU),
                expression(GB),
                expression(paste(Peers))),
       col=c("grey", "red", "blue"), 
       lty=1:1, 
       cex=0.8
)

#### ADF & KPSS tests ####
## ADF test (null = unit root)
## KPSS test (null = trend stationarity)

# ADF & KPSS - whole EU
adf.test(df_aggregate$x) # the data set has a unit root. 
kpss.test(df_aggregate$x, null = "Trend") # the data set has a unit root

# ADF & KPSS - GB
adf.test(df_GB_aggregate$x) # the data set has a unit root. The unit root null hypothesis was not rejected at alpha - 0.05
kpss.test(df_GB_aggregate$x, null = "Trend") # the data has a unit root, trend stationarity null is rejected at alpha = 0.05

# ADF & KPSS - EU excl GB
adf.test(df_exclGB_aggregate$x) # the data set has a unit root at alpha = 0.05
kpss.test(df_exclGB_aggregate$x, null = "Trend") # data has a unit at alpha = 0.10

#### Detrending ####
# detrending - EU
df_diff<- diff(df_aggregate$x)
adf.test(df_diff) # the data is stationary at alpha=0.10
kpss.test(df_diff, null = "Trend") # the data is trend stationary at alpha = 0.10
#should I use null = "Trend" or "Level"?

plot(ts(df_diff, start = 2005 , end = 2019),
     ylab = "CO2 emissions (t)", xlab = "Year",
     type = "l", 
     main = "Differenced Verified Emissions EU-ETS",
     #ylim=c(800000000, 1600000000),
     xlim=c(2005,2019))

# detrending - GB
df_diff_GB<- diff(df_GB_aggregate$x)
adf.test(df_diff_GB) # the data still has a unit root at alpha =0.05
kpss.test(df_diff_GB, null = "Trend") # doubtfull: the data is stationary at alpha = 0.05;
# the null that x is trend stationary is not rejected at alpha = 0.05, but it is rejected at alpha = 0.10 

plot(ts(df_diff_GB, start = 2005 , end = 2019),
     ylab = "CO2 emissions (t)", xlab = "Year",
     type = "l", 
     main = "Differenced Emissions in the GB, EU-ETS",
     #ylim=c(800000000, 1600000000),
     xlim=c(2005,2019))

#detrending - EU excl GB
df_diff_exclGB<- diff(df_exclGB_aggregate$x)
adf.test(df_diff_exclGB) #still has a unit root

test <- ur.df(df_diff_exclGB, type = "trend", lags=1)
summary(test)

kpss.test(df_diff_exclGB, null = "Trend") #still has a unit root at alpha = 0.10

plot(ts(df_diff_exclGB, start = 2005 , end = 2019),
     ylab = "CO2 emissions (t)", xlab = "Year",
     type = "l", 
     main = "Differenced Emissions in the EU-ETS, excl GB",
     #ylim=c(800000000, 1600000000),
     xlim=c(2005,2019))

################################################################################  
################################### 3. ARCO #################################### 
################################################################################  

##### 3.1. ArCo Linear Regressions #####
# ArCo can be fitted using linear regression functions below
# note: in the application below, we use LASSO, not this linear regression.
fn=function(X,y){
        return(lm(y~X))
}
p.fn=function(model,newdata){
        b=coef(model)
        return(cbind(1,newdata) %*% b)}

##### 3.2. Arco data prep ####
# data prep: UK (treated) and other EU countries (peers) in 1 df 
df2 <- aggregate(verified ~year+country_id, data=df, FUN=sum, na.rm=TRUE, na.action = NULL)

df2<-acast(df2, year~country_id, value.var="verified")
#df2<- subset(df2, select=-c(YEAR)) # remove year column
col_order <- c("GB","AT","BE","BG","CY","CZ","DE","DK","EE","ES","FI","FR",
               "GR","HR","HU","IE","IS","IT","LI","LT","LU","LV","MT","NL",
               "NO","PL","PT","RO","SE","SI","SK","XI")
df2 <- df2[, col_order]
df2 <- df2[-c(16:26),]

# T, N, N-1, t0
T <- nrow(df2) # T: total nr of observations, t=1,...,T
t0 <- 9 #t0: year of implementation 
#countries <- unique(df[c("country_id")])
#N <- nrow(countries) # total units i = 1,...,n
N <- ncol(df2)
n_peers <- N-1 # nr of peers

data <- list(df2) # appropriate input in ArCo fit

#### 3.3 ArCo application ####
# approximation of t0 - does not work 
t0_hat=estimate_t0(data = data, fn = fn, p.fn = p.fn, treated.unit = 1) 

# estimating ArCo  - doesn't work, probably because we do not have enough observations
# T is small compared to N-1
arco1 <- fitArCo(data = data, fn = fn, p.fn = p.fn, treated.unit = 1 , 
                t0 = 9)# VCOV.type="nw", kernel.type="TukeyHanning")
plot(arco1,plot=1,display.fitted = TRUE)
plot(arco1)

arco2 <- fitArCo(data = data, fn = cv.glmnet, p.fn = predict, treated.unit = 1 , 
                 t0 = 9)# VCOV.type="nw", kernel.type="TukeyHanning")
plot(arco2,plot=1,display.fitted = TRUE)
# does not look good either + gives a warning message

################################################################################  
################################### 4. SYNTH ################################### 
################################################################################  

#### 4.1 SC data prep ####
dfLong <- melt(df2, value.name="verified_emissions") #Var1=year, Var2=Country_code
colnames(dfLong) <- c("year", "country_code","verified_emissions")
dfLong <- transform(dfLong, unit.id=as.numeric(factor(country_code)))
dfLong$year <- as.numeric(dfLong$year) # time.variable must be a numeric variable
dfLong$country_code <- as.character(dfLong$country_code) # unit.names.variable must be a character variable
dfLong$verified_emissions <- as.numeric(dfLong$verified_emissions) # predictor must be a numeric variable

SC_dataprep.out<- dataprep(
                foo = dfLong,
                predictors = "verified_emissions",
                predictors.op = "mean", #mean is the default
                dependent = "verified_emissions",
                unit.variable = "unit.id",  # this has to be a numeric variable
                time.variable = "year",
                treatment.identifier = 1, #united kingdom
                controls.identifier = c(2:32), # you have to specify at least two control units
                time.predictors.prior = c(2005:2012),
                time.optimize.ssr = c(2005:2013), 
                unit.names.variable = "country_code", 
                time.plot = 2005:2019
                )

#### 4.2 apply SC ####
synth.out <- synth(SC_dataprep.out)

## summarizing the results by accessing the output from synth.out directly
round(synth.out$solution.w,3) # contains the unit weights
synth.out$solution.v  ## contains the predictor weights. 
# MH: here I only have 1 predictor so weight = 1

## the period by period discrepancies between the treated unit and its SC unit
gaps<- SC_dataprep.out$Y1plot-(
        SC_dataprep.out$Y0plot%*%synth.out$solution.w)

# Get result tables - V and W weights plus balance btw. treated and SC
synth.tables <- synth.tab(
        dataprep.res = SC_dataprep.out,
        synth.res = synth.out)
print(synth.tables)

## summary plots for outcome trajectories of the treated and the SC unit

#### 4.2 plot SC ####

## plot in levels (path; treated and synthetic)
path.plot(dataprep.res = SC_dataprep.out, 
          synth.res = synth.out,
          Ylab = c("Verified CO2 emissions"),
          Xlab = c("Year"), 
          #Ylim = c(0,13), 
          Legend = c("UK","Synthetic UK"),
) 
## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = SC_dataprep.out,
          synth.res = synth.out,
          Ylab = c("Gap in Verified CO2 emissions (MtCO2e)"),
          Xlab = c("Year"))
