# ArCo and SC on Production data

setwd("/cloud/project")
cat("\014")        
rm(list = ls())

##### LIBRARY LOAD ####
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

##### Load data & set variables #####
df <- read_excel("test_production_2.xlsx")
df<- as.data.frame(df)

# set certain variables
alpha = 0.05
t0=24
#treated=33

##### ARCO on GROSS ELECTRICITY PRODUCTION by SOLID FOSSIL FUELS #####
# data format for ArCo
rownames(df) <- df[,1]
data <- df[-c(1:5)]

# prepare data for ArCo
# Now, all countries are included. Countries may be removed later
data <- subset(data,select=
                 c( "United Kingdom","Austria","Belgium","Bulgaria","Croatia",
                    "Cyprus","Czechia","Denmark","Estonia","Finland",
                    "France","Germany (until 1990 former territory of the FRG)","Greece","Hungary","Ireland", "Italy",
                    "Latvia","Lithuania","Luxembourg","Malta","Netherlands", 
                    "Poland","Portugal","Romania","Slovakia","Slovenia",
                    "Spain", "Sweden", "Norway"))  
data_arco <- list(data)

arco <- fitArCo(data = data_arco,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = 1, t0 = t0, 
                boot.cf=TRUE, R=1000, l=3, 
                VCOV.type = "nw", prewhitening.kernel = TRUE)
plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Gross Electricity Production UK ArCo - Solid Fossil Fuels", ylab = "GEP (TJ)")
delta<-arco$delta 
p<-arco$p.value 
test <- p<alpha

###### UNIT ROOT & TREND STATIONARITY TESTS #####
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
tests <- unitroot_table(data) 

###### DIFFERENCING #####
#difference the data
diff_data <- diff(as.matrix(data))
diff_tests <- unitroot_table(as.data.frame(diff_data)) # only bulgaria still contains unit root
diff_data <- subset(diff_data[,-4]) #remove bulgaria

###### APPLY ARCO ON DIFF DATA #####
#apply ArCo on  differenced data
input2 <- list(diff_data)
arco2 <- fitArCo(data = input2,  fn = cv.glmnet, p.fn = predict, 
                 treated.unit = 1, t0 = 24, 
                 boot.cf=TRUE, R=1000, l=3, 
                 VCOV.type = "nw", prewhitening.kernel = TRUE)
plot(arco2, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Differenced Gross Electricity Production UK ArCo - Solid Fossil Fuels", ylab = "diff. GEP (TJ)")
delta2<-arco2$delta 
p2<-arco2$p.value 
p2<alpha #TRUE 
test2 <- (p2<alpha)

###### SC ######
# prepare data for SC
df <- df[-c(2)] # remove SIEC column
dfLong <- melt(df, id.vars="TIME", 
               variable.name = "GEO", value.name="GEP") #transform to a long df
N <- ncol(df)-1 # N = nr of countries, incl treated unit.
T <- nrow(df) # T = nr of Years
dfLong$year.id <- rep(c(1:T),times=N)
dfLong <- transform(dfLong, unit.id=as.numeric(factor(GEO)))

dfLong$TIME <- as.numeric(dfLong$TIME) # time.variable must be a numeric variable
dfLong$GEO <- as.character(dfLong$GEO) # unit.names.variable must be a character variable
dfLong$GEP <- as.numeric(dfLong$GEP) # predictor must be a numeric variable

# create matrices from panel data that provide inputs for synth()
dataprep.out<-
  dataprep(
    foo = dfLong,
    predictors = "GEP",
    predictors.op = "mean", #mean is the default
    dependent = "GEP",
    unit.variable = "unit.id",  # this has to be a numeric variable
    time.variable = "TIME",
    # not sure what the below command is for.
    # special.predictors = list( 
    #  list("Y", 1991, "mean"),
    # list("Y", 1985, "mean"),
    # list("Y", 1980, "mean")
    # ),
    treatment.identifier = 33, #united kingdom
    controls.identifier = c(4:32), # you have to specify at least two control units
    time.predictors.prior = c(1990:2012),
    time.optimize.ssr = c(1990:2013), 
    unit.names.variable = "GEO", 
    time.plot = 1990:2018
  )

# run synth
# identifies weights that create the best possible synthetic control unit for the treated.
synth.out <- synth(dataprep.out)

## summarizing the results by accessing the output from synth.out directly
round(synth.out$solution.w,3) # contains the unit weights
synth.out$solution.v  ## contains the predictor weights. 
# MH: here I only have 1 predictor so weight = 1

gaps<- dataprep.out$Y1plot-(
  dataprep.out$Y0plot%*%synth.out$solution.w
) ; gaps

# Get result tables - V and W weights plus balance btw. treated and SC
synth.tables <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res = synth.out)
print(synth.tables)

## summary plots for outcome trajectories of the treated and the SC unit

## plot in levels (path; treated and synthetic)
path.plot(dataprep.res = dataprep.out, 
          synth.res = synth.out,
          Ylab = c("GEP (TJ)"),
          Xlab = c("Year"), 
          #Ylim = c(0,13), 
          Legend = c("UK","Synthetic UK"),
          Main = "Gross Electricity Production UK SC - Solid Fossil Fuels"
) 
## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep.out,
          synth.res = synth.out,
          Ylab = c("Gap in GEP (TJ)"),
          Xlab = c("Year"),
          Main = "Gaps between GEP in UK and its SC - Solid Fossil Fuels")