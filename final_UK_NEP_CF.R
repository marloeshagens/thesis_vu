################################################################################
# FINAL VERSION APPLICATION CODE 
# July 16 2021
# author: Marloes Hagens
################################################################################
# application of counterfactual methods ArCo and SC
# for the case of a carbon tax with the aim to decrease CO2 emissions 
# in the electricity sector of the treated unit the United Kingdom
# and EU countries as donor pool countries
# using net electricity production (NEP) data from Eurostat
################################################################################
# table of contents #
# 1 Load Libraries and Functions
# 2 Load and Prepare Data
# 3 Unit Root Tests
# 4 Set Variables
# 5 Counterfactuals ArCo and SC
# 6 Additional Descriptive Plots
################################################################################

#### 1 Load Libraries and Functions ####
setwd("/cloud/project")
# setwd("~/thesis_vu-master")
cat("\014")        
rm(list = ls())

library(tidyr) # for "spread" function in wide df
library(ArCo)
library(tseries)
library(Synth)
library(glmnet) # cv.glmnet (LASSO)
library(reshape) # function melt()
library(directlabels)
library(ggplot2)

set.seed(42)

ptm <- proc.time() # start timer

#### 2 Load and Prepare Data ####
# read data  
data <- read.csv("nrg_ind_peh__custom_1070025_20210616_221539.sdmx.csv")
df <- subset(data, select=-c(DATAFLOW,freq, unit, OBS_FLAG))
df <- subset(df, (nrg_bal == "NEP"  # select only Net Electricity Production
                  & plants == "ELC" # select only electricity production units
                  & operator == "PRR_MAIN" # select only operators which main 
                  # activity is electricity production
))
df <- subset(df, select=-c(nrg_bal,plants, operator))
# exclude aggregate units
df2 <- subset(df, (geo != "EA19" & geo!= "EU27_2020" & geo!="EU28")) 
dfwide<-spread(df2, geo, OBS_VALUE)
CF <- subset(dfwide, siec == "CF", select=-siec) 
# disregard columns that have NA values
CF <- CF[,colSums(is.na(CF)) == 0] 
# disregard countries with 0 values
CF <- subset(CF, select=-c(EE,LT,LV,PL,RS,UA)) 

# define T, N
T <- nrow(CF)
N <- ncol(CF) - 1

# ArCo data format - non-differenced and differenced
data_arco <- CF
rownames(data_arco) <- data_arco[,1]
data_arco <- data_arco[-c(1)]
df_arco <- list(data_arco) # non-differenced data for the ArCo 
diff_data <- diff(as.matrix(data_arco))
diff_df_arco<-list(diff_data) # differenced data for the ArCo input

# SC data format - non-differenced and differenced
# transform data to a long df
df_sc <- melt(CF, id.vars="TIME_PERIOD", variable.name = "geo", value.name="NEP") 
colnames(df_sc) <- c("TIME_PERIOD", "geo", "NEP") # correct column names
df_sc$year.id <- rep(c(1:T),times=N)

# non-differenced data for the SC 
df_sc <- transform(df_sc, unit.id=as.numeric(factor(geo))) 
df_sc$geo <- as.character(df_sc$geo) # must be character 
df_sc$NEP <- as.numeric(df_sc$NEP) # must be numeric 
df_sc$TIME_PERIOD <- as.numeric(df_sc$TIME_PERIOD) # must be numeric 

# differenced data for the SC
diff_dfwide <- diff_data
diff_dfwide <- as.data.frame(cbind(c(1991:2019), diff_dfwide))
colnames(diff_dfwide) <- colnames(CF)
# transform data to a long df
diff_df_sc <- melt(diff_dfwide, id.vars="TIME_PERIOD", variable.name = "geo", 
                   value.name="NEP")
colnames(diff_df_sc) <- c("TIME_PERIOD", "geo", "NEP") # correct column names
diff_df_sc$year.id <- rep(c(2:T),times=N) # start at 2 here, bcs of 1st diff
diff_df_sc <- transform(diff_df_sc, unit.id=as.numeric(factor(geo)))
diff_df_sc$TIME_PERIOD <- as.numeric(diff_df_sc$TIME_PERIOD) # must be numeric 
diff_df_sc$geo <- as.character(diff_df_sc$geo) # must be character 
diff_df_sc$NEP <- as.numeric(diff_df_sc$NEP) # must be numeric 

#### 3 Unit Root Tests ####
# unit root, trend and level stationarity tests
adf <- t(sapply(data_arco, adf.test))
kpss.trend <- t(sapply(data_arco, kpss.test, null = "Trend"))
kpss.level <- t(sapply(data_arco, kpss.test, null = "Level"))

#### 4 Set Variables #### 

# select row that contains 2013 (= treatment year)
t0 = which(rownames(data_arco)=="2013")
t0_diff = t0-1

# select column that contains UK (= treated unit)
treated = which(colnames(data_arco)=="UK") 

# amount of bootstrap steps in the ArCo to construct confidence intervals
R = 10000 

# define donor pool
donorpool <- rep(c(1:N)) 
donorpool <- donorpool[c(-treated)]

#### 5 Counterfactuals ArCo and SC ####

# ArCo 
# ArCo - normal data
arco <- fitArCo(data = df_arco,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = treated, t0 = t0, 
                boot.cf=TRUE, R=R,
                VCOV.type = "iid")
plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Net Electricity Production UK ArCo - Combustible Fuels (VCOV type = iid)", ylab = "NEP (GWH)")
delta <- arco$delta 
p <- arco$p.value 
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

arco_nw <- fitArCo(data = df_arco,  fn = cv.glmnet, p.fn = predict, 
                   treated.unit = treated, t0 = t0, 
                   boot.cf=TRUE, R=R, 
                   VCOV.type = "nw" 
)
plot(arco_nw, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Net Electricity Production UK ArCo - Combustible Fuels (VCOV type = NW)", ylab = "NEP (GWH)")
delta_nw <- arco_nw$delta 
p_nw <- arco_nw$p.value 
test_nw <- p_nw < 0.05

# ArCo - first differenced data
t0_diff = t0-1
arco_diff <- fitArCo(data = diff_df_arco,  fn = cv.glmnet, p.fn = predict, 
                     treated.unit = treated, t0 = t0_diff, 
                     boot.cf=TRUE, R=R, 
                     VCOV.type = "iid"
)
plot(arco_diff, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "ΔNEP UK ArCo - Combustible Fuels (VCOV type = iid)", ylab = "ΔNEP (GWH)")
delta_diff <- arco_diff$delta 
p_diff <- arco_diff$p.value 
test_diff <- p_diff<0.05

e_diff<-arco_diff$residuals 
par(mfrow=c(2,2), mar=c(2,3.9,2,3.9))
plot(ts(e_diff, 
        frequency = 1), type = "l", cex = 0.5, 
     ylab = "residuals")
hist(e_diff, prob = TRUE, main = "")
lines(density(e_diff))
qqnorm(e_diff, main = "")
qqline(e_diff)
acf(e_diff, lag.max = 12, main = "")

arco_diff_nw <- fitArCo(data = diff_df_arco,  fn = cv.glmnet, p.fn = predict, 
                        treated.unit = treated, t0 = t0_diff, 
                        boot.cf=TRUE, R=R, 
                        VCOV.type = "nw"
)
plot(arco_diff_nw, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "ΔNEP UK ArCo - Combustible Fuels (VCOV type = NW)", ylab = "ΔNEP (GWH)")
delta_diff_nw <- arco_diff_nw$delta 
p_diff_nw <- arco_diff_nw$p.value 
test_diff_nw <- p_diff_nw<0.05

# SC
# SC - normal data
dataprep.out<- dataprep(foo = df_sc,# create matrices from panel data that provide inputs for synth()
                        predictors = "NEP",
                        predictors.op = "mean", # mean is the default
                        dependent = "NEP",
                        unit.variable = "unit.id",  # this has to be a numeric variable
                        time.variable = "TIME_PERIOD",
                        treatment.identifier = treated, 
                        controls.identifier = donorpool,
                        time.predictors.prior = c(1990:2012), 
                        time.optimize.ssr = c(1990:2013), 
                        unit.names.variable = "geo",
                        time.plot = c(1990:2019))

# run synth
synth.out <- synth(dataprep.out) # identifies weights that create the best possible synthetic control unit for the treated.
round(synth.out$solution.w,3) # # SC results summary, contains the unit weights # 
synth.out$solution.v  # contains the predictor weights: I only use 1 predictor so weight = 1 #

synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
print(synth.tables)

# SC - first differenced data
diff.dataprep.out<- dataprep(
  foo = diff_df_sc,
  predictors = "NEP",
  predictors.op = "mean", # mean is the default
  dependent = "NEP",
  unit.variable = "unit.id",  # this has to be a numeric variable
  time.variable = "TIME_PERIOD",
  treatment.identifier = treated, 
  controls.identifier = donorpool,
  time.predictors.prior = c(1991:2012),
  time.optimize.ssr = c(1991:2013), 
  unit.names.variable = "geo", 
  time.plot = c(1991:2019)
)

# run synth
diff.synth.out <- synth(diff.dataprep.out) # identifies weights that create the best possible synthetic control unit for the treated.
round(diff.synth.out$solution.w,3) # contains the unit weights # 
diff.synth.out$solution.v  # contains the predictor weights: I only use 1 predictor so weight = 1 #

diff.synth.tables <- synth.tab(dataprep.res = diff.dataprep.out, synth.res = diff.synth.out)
print(diff.synth.tables)

# Counterfactuals
Y <- dataprep.out$Y1plot # original data 
Y_diff<-diff.dataprep.out$Y1plot# original first diff. data

arco_cf <- c(fitted(arco_nw)[,1], arco_nw$cf[,1])
sc_cf <- dataprep.out$Y0plot%*%synth.out$solution.w

arco_cf_diff <- c(fitted(arco_diff_nw)[,1], arco_diff_nw$cf[,1])
sc_cf_diff <- diff.dataprep.out$Y0plot%*%diff.synth.out$solution.w

par(mfrow=c(1,1))
plot(Y,type="l",ylab="NEP (GWH)",xlab="Time (Years)",
     main = "NEP (combustible fuels) for the UK and its counterfactuals")
lines(sc_cf,col=2)
abline(v=24,col=8,lty=2)
lines(arco_cf, col=4)
legend("bottomleft",legend=c("UK","SC UK","ArCo UK"),
       col=c(1,2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

plot(Y_diff,type="l",ylab="ΔNEP (GWH)",xlab="Time (Years)", 
     main = "ΔNEP (combustible fuels) for the UK and its counterfactuals") # (CF) = combustible fuels
lines(sc_cf_diff,col=2)
abline(v=23,col=8,lty=2)
lines(arco_cf_diff,col=4)
legend("bottomleft",legend=c("UK","SC UK","ArCo UK"),
       col=c(1,2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

# calculate the gap between Y and its estimated counterfactual values 
gap_arco <- Y - arco_cf
gap_sc <- Y - sc_cf
gap_arco_diff <-  Y_diff - arco_cf_diff
gap_sc_diff <- Y_diff - sc_cf_diff

plot(gap_arco,type="l",ylab="NEP - δ (GWH)",xlab="Time (Years)", col = 4, ylim=c(-80000,40000), 
     main = "Gaps between the NEP and its counterfactuals (δ) in UK")
lines(gap_sc,col=2)
abline(h=0)
abline(v=24,col=8,lty=2)
legend("bottomleft",legend=c("NEP - δ_t Synth","NEP - δ_t ArCo"),
       col=c(2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

plot(gap_arco_diff,type="l",ylab="ΔNEP - δ_t (GWH)",xlab="Time (Years)", 
     col = 4, ylim=c(-30000,20000), 
     main = "Gaps between the ΔNEP and its counterfactuals (δ) in UK")
lines(gap_sc_diff,col=2)
abline(h=0)
abline(v=24,col=8,lty=2)
legend("bottomleft",legend=c("ΔNEP - δ_t Synth","ΔNEP - δ_t ArCo"),
       col=c(2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

# compare ArCo and SC magnitude of values (incl shapiro-wilk, t-test)
par(mfrow=c(1,1), mar=c(2,3.9,2,3.9))
d <- ts((arco_cf - sc_cf), start = c(1990,1), end = 2019)
plot(d, type="l", main=expression("ΔArCo-SC"), ylab=expression("ΔY"), xlab="Time (Years)")
abline(v=2013,col=8,lty=2)
abline(h=0,lty=2)

# Shapiro-Wilk normality tests #
shapiro.test(d)
shapiro.test(arco_cf)
shapiro.test(sc_cf)
# normality cannot be assumed - no t test possible

d_diff <- ts((arco_cf_diff - sc_cf_diff), start = c(1990,2), end = 2019)
plot(d_diff,type="l", main=expression("ΔArCo-SC (first differenced data)"), ylab=expression("ΔY"), xlab="Time (Years)")
abline(v=2013,col=8,lty=2)
abline(h=0,lty=2)

shapiro.test(d_diff) 
shapiro.test(arco_cf_diff) 
shapiro.test(sc_cf_diff) 
# normality can be assumed

t.test(arco_cf_diff, sc_cf_diff) # not significant
# t.test(d_diff, mu=0) # not significant

# Sum of gaps after T0
# thus: total amount of electricity (GWH) that would have been produced extra, 
# according to the ArCo and SC
pre_t0_gaps <- cbind(sum((gap_arco[1:23,])), (sum(gap_sc[1:23,])), 
                     sum((gap_arco_diff[1:22,])), sum((gap_sc_diff[1:22,])))
colnames(pre_t0_gaps)<-c("gap_arco","gap_sc","gap_arco_diff","gap_sc_diff")
rownames(pre_t0_gaps)<- "sum"

post_t0_gaps <- cbind(sum((gap_arco[24:30,])), (sum(gap_sc[24:30,])), 
                      sum((gap_arco_diff[23:29,])), sum((gap_sc_diff[23:29,])))
colnames(post_t0_gaps)<-c("gap_arco","gap_sc","gap_arco_diff","gap_sc_diff")
rownames(post_t0_gaps)<- "sum"

# ATET
DeltaT<- cbind(mean((gap_arco[24:30,])), (mean(gap_sc[24:30,])))
colnames(DeltaT)<-c("Delta_T ArCo", "Delta_T SC")
DeltaT_diff<- cbind(mean((gap_arco_diff[23:29,])), mean((gap_sc_diff[23:29,])))
colnames(DeltaT_diff)<-c("Delta_T (diff) ArCo", "Delta_T (diff) SC")    

#### 6 Additional Descriptive Plots ####
# exclude siec's with very low values
df_EU <- subset(df, ( geo=="EU28" &   siec !="X9900H" & siec !="X9900" & siec !="RA500" & siec !="RA130" & siec !="RA200"))
df_GB <- subset(df, ( geo=="UK" &   siec !="X9900H" & siec !="X9900" & siec !="RA500" & siec !="RA130" & siec !="RA200"))

EU_siec_plot <- ggplot(df_EU, aes(x = TIME_PERIOD, y = OBS_VALUE, group = siec, 
                                  colour = siec)) + 
  geom_line(aes(colour = siec)) 
print(EU_siec_plot 
      + labs(title= "Net Electricity Production in EU28 countries",
             y="Net Electricity Production (GWH)", x = "Time (Years)")
      + theme_light())

UK_siec_plot <- ggplot(df_GB, aes(x = TIME_PERIOD, y = OBS_VALUE, group = siec, colour = siec)) + 
  geom_line(aes(colour = siec)) 
print(UK_siec_plot 
      + labs(title= "Net Electricity Production in the United Kingdom",
             y="Net Electricity Production (GWH)", x = "Time (Years)")
      + theme_light())

# for combustible fuels use only, 
# select a few (arbitrary) countries to show trends
df_CF_UK <- subset(df, ( siec=="CF" & geo=="UK"))
df_CF_EU28 <- subset(df, ( siec=="CF" & geo=="EU28"))
df_CF_DE <- subset(df, ( siec=="CF" & geo=="DE"))
df_CF_FR <- subset(df, ( siec=="CF" & geo=="FR"))
df_CF_IE <- subset(df, ( siec=="CF" & geo=="IE"))

# normalizing function, such we can compare the data to each
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) 
}
df_CF_UK$ynorm <- normalize(df_CF_UK$OBS_VALUE)
df_CF_EU28$ynorm <- normalize(df_CF_EU28$OBS_VALUE)
df_CF_DE$ynorm <- normalize(df_CF_DE$OBS_VALUE)
df_CF_FR$ynorm <- normalize(df_CF_FR$OBS_VALUE)
df_CF_IE$ynorm <- normalize(df_CF_IE$OBS_VALUE)
normplotdata <- rbind(df_CF_UK, df_CF_EU28, df_CF_DE, 
                      df_CF_FR, df_CF_IE)
colnames(normplotdata) <- c("siec", "Country","Year","NEP","Norm.value")

normplot<- ggplot(normplotdata, aes(x = Year, y = Norm.value, group = Country))+
  geom_line(aes(linetype = Country, color = Country))+
  geom_point(aes(color = Country))+
  theme(legend.position = "top")
normplot + theme_light() + labs(y= "Normalized NEP (GWH)", x = "Year") 

NEPplot<- ggplot(normplotdata, aes(x = Year, y = NEP, group = Country))+
  geom_line(aes(linetype = Country, color = Country))+
  geom_point(aes(color = Country))+
  theme(legend.position = "top")
NEPplot + theme_light() + labs(y= "NEP (GWH)", x = "Year") 

timer <- proc.time() - ptm # stop timer
