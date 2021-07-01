setwd("/cloud/project")
# setwd("~/thesis_vu-master")
cat("\014")        
rm(list = ls())

library(tidyr) # for "spread" function in wide df
library(ArCo)
library(tseries)
library(Synth)

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
CF <- subset(dfwide, siec == "CF", select=-siec) # only use CF for now.

# N9000<- subset(dfwide, siec == "N9000", select=-siec)
# RA100 <- subset(dfwide, siec == "RA100", select=-siec)
# RA300 <- subset(dfwide, siec == "RA300", select=-siec)
# RA400 <- subset(dfwide, siec == "RA400", select=-siec)
# RA500 <- subset(dfwide, siec == "RA500", select=-siec)
# X9900 <- subset(dfwide, siec == "X9900", select=-siec)
# X9900H <- subset(dfwide, siec == "X9900H", select=-siec)
# TOTAL <- subset(dfwide, siec == "TOTAL", select=-siec)
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

# data cleaning
CF <- CF[,colSums(is.na(CF)) == 0] # disregards columns that have NA values
CF <- subset(CF, select=-c(EE,LT,LV,PL,RS,UA)) #disregard countries with 0 values

#further cleaning: discuss with richard first:
# CF <- subset(CF, select=-c(SE, NO, AL, IS)) #disregard countries with too low CF use (from beginning onwards)

# ArCo compatible data 
df_arco <- CF
rownames(df_arco) <- df_arco[,1]
df_arco <- df_arco[-c(1)]
data_arco <- list(df_arco)

# define T, N, t0, treated unit
T <- nrow(CF)
N <- ncol(CF) - 1
t0 = which(rownames(df_arco)=="2013") # year 2013 
tprior <- t0-1
treated = which(colnames(df_arco)=="UK") # select column that contains UK (= treated unit)

# ARCO on UK's Net Electricity Production for Combustible Fuels in Gigawatt-hour #
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
# moreover, the delta's are in each other LB and UB. Therefore, we keep the iid situation

# UNIT ROOT & TREND STATIONARITY TESTS #
adf <- lapply(df_arco, adf.test)
kpss <- lapply(df_arco, kpss.test, null = "Trend")
unitroot_table<- function(df){
  p <- ncol(df)
  df_stats <- data.frame(var=names(df),
                         adf.pvalue=sapply(df, function(v) adf.test(ts(v), alternative = "stationary")$p.value),
                         kpss.pvalue=sapply(df, function(v) kpss.test(ts(v), null = "Trend")$p.value)
  )
  df_stats$unitroot_adf <- ifelse(df_stats$adf.pvalue < 0.05, 0, 1)
  df_stats$unitroot_kpss <- ifelse(df_stats$kpss.pvalue < 0.05, 1, 0)
  df_stats$conclusion <- ifelse(df_stats$unitroot_kpss + df_stats$unitroot_adf > 1, 
                                "contains unit root", "trend stationary or unclear")
  row.names(df_stats) <- c()
  df_stats
}
tests <- unitroot_table(df_arco) 

# difference the data #
# diff_data <- diff(as.matrix(log(df_arco)))
diff_data <- diff(as.matrix(df_arco))
list_diff_data<-list(diff_data)

# arco on diff data #
t0_diff = t0-1
arco_diff <- fitArCo(data = list_diff_data,  fn = cv.glmnet, p.fn = predict, 
                 treated.unit = treated, t0 = t0_diff, 
                 boot.cf=TRUE, R=500, 
                 VCOV.type = "nw"
                 #prewhitening.kernel = TRUE
                 )
plot(arco_diff, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main= "Differenced NEP UK ArCo - Combustible Fuels (VCOV type = NW)", ylab = "differenced NEP (GWH)")
delta_diff<-arco_diff$delta 
p_diff<-arco_diff$p.value 
test_diff <- p_diff<0.05
# the differenced data is compared to the nw version of the normal data arco
# both the iid and nw version give a significant result
# only if I use prewhitning kernel = TRUE at differenced data, it gives a non-negative result.
# BUT WHY

# SC #

# prepare data for SC
df_sc <- melt(CF, id.vars="TIME_PERIOD", 
                variable.name = "geo", value.name="NEP") #transform to a long df
df_sc$year.id <- rep(c(1:T),times=N)
df_sc <- transform(df_sc, unit.id=as.numeric(factor(geo)))

df_sc$TIME_PERIOD <- as.numeric(df_sc$TIME_PERIOD) # time.variable must be a numeric variable
df_sc$geo <- as.character(df_sc$geo) # unit.names.variable must be a character variable
df_sc$NEP <- as.numeric(df_sc$NEP) # predictor must be a numeric variable

# define donorpool #
donorpool <- rep(c(1:N))
donorpool <- donorpool[c(-treated)]

# define time periods
time_prior = c(1990:2012)
time_prior_and_t0 = c(1990:2013)
total_time = c(1990:2019)

# create matrices from panel data that provide inputs for synth()
dataprep.out<- dataprep(foo = df_sc,
predictors = "NEP",
predictors.op = "mean", # mean is the default
dependent = "NEP",
unit.variable = "unit.id",  # this has to be a numeric variabletime.variable = "TIME_PERIOD",
treatment.identifier = treated, 
controls.identifier = donorpool,
time.predictors.prior = time_prior, 
time.optimize.ssr = time_prior_and_t0, 
unit.names.variable = "geo",
time.plot = total_time)

# run synth
synth.out <- synth(dataprep.out) # identifies weights that create the best possible synthetic control unit for the treated.
# X0 the control cases after the treatment
# X1 the control case before the treatment
# Z1 the treatment case before the treatment
# Z0 the treatment case after the treatment

# SC results summary #
round(synth.out$solution.w,3) # contains the unit weights # 
synth.out$solution.v  # contains the predictor weights: I only use 1 predictor so weight = 1 #

# calculate difference in trend between data and its synthetic control #
# gaps<- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%synth.out$solution.w)

# pre built tables from synth objects #
synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
print(synth.tables)

# summary plots for outcome trajectories of the treated and the SC unit #
# plot treatment vs control outcomes for pre and post periods # 
# dev.off() 
# path.plot(dataprep.res = dataprep.out, 
#           synth.res = synth.out,
#           Ylab = c("NEP (GWH)"),
#           Xlab = c("Year"), 
#           #Ylim = c(0,13), 
#           Legend = c("UK","Synthetic UK"),
#           Main = "Net Electricity Production UK SC - Combustible Fuels"
# ) 
# # plot the gaps (treated - synthetic) # 
# gaps.plot(dataprep.res = dataprep.out,
#           synth.res = synth.out,
#           Ylab = c("Gap in NEP (GWH)"),
#           Xlab = c("Year"),
#           Main = "Gaps between NEP in UK and its SC - Combustible Fuels")

### SCdiff ### 
diff_dfwide <- diff_data
diff_dfwide <- as.data.frame(cbind(c(1991:2019), diff_dfwide))
colnames(diff_dfwide) <- colnames(CF)

diff_df_sc <- melt(diff_dfwide, id.vars="TIME_PERIOD", 
              variable.name = "geo", value.name="NEP") #transform to a long df
diff_df_sc$year.id <- rep(c(2:T),times=N) # note we start at 2 here, because of the differenced data
diff_df_sc <- transform(diff_df_sc, unit.id=as.numeric(factor(geo)))

diff_df_sc$TIME_PERIOD <- as.numeric(diff_df_sc$TIME_PERIOD) # time.variable must be a numeric variable
diff_df_sc$geo <- as.character(diff_df_sc$geo) # unit.names.variable must be a character variable
diff_df_sc$NEP <- as.numeric(diff_df_sc$NEP) # predictor must be a numeric variable

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
# diff.gaps<- diff.dataprep.out$Y1plot-( # calculate difference in trend between data and its synthetic control #
#   diff.dataprep.out$Y0plot%*%diff.synth.out$solution.w
# )
# diff.gaps

# pre built tables from synth objects #
diff.synth.tables <- synth.tab(dataprep.res = diff.dataprep.out, synth.res = diff.synth.out)
print(diff.synth.tables)

# plotted together with arco below
# summary plots for outcome trajectories of the treated and the SC unit #
# plot treatment vs control outcomes for pre and post periods # 
# path.plot(dataprep.res = diff.dataprep.out, synth.res = diff.synth.out,
#           Ylab = c("Diff NEP (GWH)"), Xlab = c("Year"),Legend = c("UK","Synthetic UK"),
#           Main = "Net Electricity Production UK SC - Combustible Fuels"
# ) 

# # plot the gaps (treated - synthetic) # 
# gaps.plot(dataprep.res = diff.dataprep.out, synth.res = diff.synth.out,
#           Ylab = c("Gap in differenced NEP (GWH)"), Xlab = c("Year"),
#           Main = "Gaps between differenced NEP in UK and its SC - Combustible Fuels")

# dev.off()
# NORMAL DATA SC AND ARCO in 1 plot
arco_cf <- c(fitted(arco_)[,1], arco_$cf[,1])
sc_cf <- dataprep.out$Y0plot%*%synth.out$solution.w
Y <- dataprep.out$Y1plot

plot(Y,type="l",ylab="Y (GWH)",xlab="Time (Years)",
     main = "NEP (CF) for the UK and the UK's counterfactuals") # (CF) = combustible fuels
lines(sc_cf,col=2)
abline(v=24,col=8,lty=2)
lines(arco_cf, col=4)
legend("bottomleft",legend=c("UK","SC UK","ArCo UK"),
       col=c(1,2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

# DIFFERENCED DATA SC AND ARCO in 1 plot
arco_cf_diff <- c(fitted(arco_diff)[,1], arco_diff$cf[,1])
sc_cf_diff <- diff.dataprep.out$Y0plot%*%diff.synth.out$solution.w
Y_diff<-diff.dataprep.out$Y1plot

plot(Y_diff,type="l",ylab="Y (GWH)",xlab="Time (Years)", 
     main = "Differenced NEP (CF) for the UK and the UK's counterfactuals") # (CF) = combustible fuels
lines(diff_sc_cf,col=2)
abline(v=23,col=8,lty=2)
lines(diff_arco_cf,col=4)
legend("bottomleft",legend=c("UK","SC UK","ArCo UK"),
       col=c(1,2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

# calculate the gap between estimated counterfactual magnitudes/values and Y
gap_arco <- arco_cf-Y
gap_sc <- sc_cf-Y
plot(gap_arco,type="l",ylab="CF - Y (GWH)",xlab="Time (Years)", col = 4, ylim=c(-40000,80000), 
     main = "Gaps between counterfactuals and the NEP in UK")
lines(gap_sc,col=2)
abline(h=0)
abline(v=24,col=8,lty=2)
legend("bottomleft",legend=c("Y - Synth","Y - ArCo"),
       col=c(2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

gap_arco_diff <- arco_cf_diff - Y_diff
gap_sc_diff <- sc_cf_diff - Y_diff
plot(gap_arco_diff,type="l",ylab="CF - ΔY (GWH)",xlab="Time (Years)", col = 4, ylim=c(-20000,30000), 
     main = "Gaps between counterfactuals and the first difference NEP in UK")
lines(gap_sc_diff,col=2)
abline(h=0)
abline(v=24,col=8,lty=2)
legend("bottomleft",legend=c("Y diff - Synth","Y diff - ArCo"),
       col=c(2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

# calculate sum of gaps after t0 
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

# COMPARE SC AND ARCO VALUES # 
# difference between the pairs
d <- ts((arco_cf - sc_cf), start = c(1990,1), end = 2019)
plot(d, type="l", main=expression("ΔArCo-SC"), ylab=expression("ΔY"), xlab="Time (Years)")
abline(v=2013,col=8,lty=2)
abline(h=0,lty=2)

d_diff <- ts((diff_arco_cf - diff_sc_cf), start = c(1990,2), end = 2019)
plot(d_diff,type="l", main=expression("ΔArCo-SC (first differenced data)"), ylab=expression("ΔY"), xlab="Time (Years)")
abline(v=2013,col=8,lty=2)
abline(h=0,lty=2)

#######
#REMOVE BELOW!
# Shapiro-Wilk normality tests #
# shapiro.test(d) 
#  the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. 
# In other words, we cannot assume the normality.
# and thus we cannot perform a statistical test comparing the two time series, whether they are significantly different from 0

# shapiro.test(d_diff) 
# W = 0.93229, p-value = 0.06308
# the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
# we can assume the normality and thus perform a statistical test comparing the two time series, whether they are significantly different from 0

# mu <- mean(d_diff)
# stdev <- sqrt(var(d_diff))