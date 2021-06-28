setwd("/cloud/project")
cat("\014")        
rm(list = ls())

library(tidyr) # for "spread" function in wide df
library(ArCo)

set.seed(42)

#### load and prepare data ####
data <- read.csv("nrg_ind_peh__custom_1070025_20210616_221539.sdmx.csv")
df <- subset(data, select=-c(DATAFLOW,freq, unit, OBS_FLAG))
df <- subset(df, (nrg_bal == "NEP"  # select only Net Electricity Production
                  & plants == "ELC" # select only electricity production units
                  & operator == "PRR_MAIN" # select only operators which main activity is electricity production
                  ))
df <- subset(df, select=-c(nrg_bal,plants, operator))
df <- subset(df, (geo != "EA19" & geo!= "EU27_2020" & geo!="EU28"))

# unique(df$geo)

#ArCo data format
dfwide<-spread(df, geo, OBS_VALUE)
CF <- subset(dfwide, siec == "CF", select=-siec)
# N9000<- subset(dfwide, siec == "N9000", select=-siec)
# RA100 <- subset(dfwide, siec == "RA100", select=-siec)
# RA300 <- subset(dfwide, siec == "RA300", select=-siec)
# RA400 <- subset(dfwide, siec == "RA400", select=-siec)
# RA500 <- subset(dfwide, siec == "RA500", select=-siec)
# X9900 <- subset(dfwide, siec == "X9900", select=-siec)
# X9900H <- subset(dfwide, siec == "X9900H", select=-siec)
# TOTAL <- subset(dfwide, siec == "TOTAL", select=-siec)
rownames(CF) <- CF[,1]
CF <- CF[-c(1)]
# rownames(N9000) <- N9000[,1]
# N9000 <- N9000[-c(1)]
# rownames(RA100) <- RA100[,1]
# RA100 <- RA100[-c(1)]
# rownames(RA300) <- RA300[,1]
# RA300 <- RA300[-c(1)]
# rownames(RA400) <- RA400[,1]
# RA400 <- RA400[-c(1)]
# rownames(RA500) <- RA500[,1]
# RA500 <- RA500[-c(1)]
# rownames(X9900) <- X9900[,1]
# X9900 <- X9900[-c(1)]
# rownames(X9900H) <- X9900H[,1]
# X9900H <- X9900H[-c(1)]

# data_arco <- list(CF, N9000, RA100, RA300, RA400, RA500, X9900, X9900H)
data<-CF[,colSums(is.na(CF)) == 0] # disregards columns that have NA values
data <- subset(data, select=-c(EE,LT,LV,PL,RS,UA)) #disregard countries with 0 values
T <- nrow(CF)
N <- ncol(CF)
t0 = which(rownames(data)=="2013") # year 2013 
treated = which(colnames(data)=="UK")

data_arco <- list(data)


# differenced data
diffdata=leveldata
diffdata$gdpcap=diff(log(diffdata$gdpcap),1)
diffdata$gdpcap=rbind(NA,diffdata$gdpcap)


#### ARCO on UK's Net Electricity Production for Combustible Fuels in Gigawatt-hour ####
arco <- fitArCo(data = data_arco,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = treated, t0 = t0, 
                boot.cf=TRUE, R=500,
                VCOV.type = "iid")
plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Net Electricity Production UK ArCo - Combustible Fuels", ylab = "NEP (GWH)")
delta<-arco$delta 
p<-arco$p.value 
test <- p<0.05
# Diagnostic Plots for the Residuals # 
  e<-arco$residuals 
  par(mfrow=c(2,2), mar=c(2,3.9,2,3.9))
  plot(ts(e, 
          frequency = 1), type = "l", cex = 0.5, 
       ylab = "residuals")
  hist(e, prob = TRUE, main = "")
  lines(density(e))
  qqnorm(e, main = "")
  qqline(e)
  acf(e, lag.max = 12, main = "")

# new arco 
arco_ <- fitArCo(data = data_arco,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = treated, t0 = t0, 
                boot.cf=TRUE, R=500, 
                VCOV.type = "nw", 
                prewhitening.kernel = TRUE)
plot(arco_, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Net Electricity Production UK ArCo - Combustible Fuels (VCOV type = NW)", ylab = "NEP (GWH)")
delta_<-arco_$delta 
p_<-arco_$p.value 
test_ <- p<0.05
# doesn't change anything as cmpared to the iid situation - so fine to keep it like this. 


#unit root tests 
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

diff_data <- diff(as.matrix(data))
list_diff_data<-list(diff_data)
# arco on diff data
t0_diff = t0-1
arco_diff <- fitArCo(data = list_diff_data,  fn = cv.glmnet, p.fn = predict, 
                 treated.unit = treated, t0 = t0_diff, 
                 boot.cf=TRUE, R=500, 
                 VCOV.type = "nw", 
                 prewhitening.kernel = TRUE)
plot(arco_diff, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Differenced NEP UK ArCo - Combustible Fuels (VCOV type = NW)", ylab = "differenced NEP (GWH)")
delta_diff<-arco_diff$delta 
p_diff<-arco_diff$p.value 
test_diff <- p<0.05
# doesn't change anything as cmpared to the iid situation - so fine to keep it like this.


###### SC ######
# prepare data for SC
#df <- df[-c(2)] # remove SIEC column
dfLong <- melt(CF, id.vars="TIME", 
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
    time.plot = 1990:2019
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
          Ylab = c("NEP (GWH)"),
          Xlab = c("Year"), 
          #Ylim = c(0,13), 
          Legend = c("UK","Synthetic UK"),
          Main = "Net Electricity Production UK SC - Combustible Fuels"
) 
## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep.out,
          synth.res = synth.out,
          Ylab = c("Gap in NEP (GWH)"),
          Xlab = c("Year"),
          Main = "Gaps between NEP in UK and its SC - Combustible Fuels")
