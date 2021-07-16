##############################################################################
### ###  NONSTATIONARY TS COUNTERFACTUAL  -  MC SIMULATION EXPERIMENT  ### ###
##############################################################################
# July 6 2021 
# author: M Hagens

setwd("~/Documents/thesis_vu-master")
cat("\014")        
rm(list = ls())

##############################################################################
############# 0. install & load packages, load functions & Rdata ############# 
##############################################################################

# install.packages(ArCo) 
# install.packages(Synth) 
# install.packages(glmnet) 
# install.packages(reshape2) 
# install.packages(Matrix) 
# install.packages(xtable)
# install.packages(writexl)

library(ArCo) # ArCo package #
library(Synth) # Synthetic Control package #
library(glmnet) # contains LASSO #
library(reshape2) # contains melt(), when using SC #
library(Matrix) # used for colSums #
library(xtable)
library(writexl)

# source self-made functions
source("randomwalk.R") # data generating process function # 
source("nonstationary_ts_counterfactual.R") # function for the nonstationary time series counterfactual experiment #
source("percentagesignificant.R") # function that computes the percentages of arco outputs that are significant
source("arcogaps.R") # function that extracts data from arco output and construct gaps 

set.seed(42)

# load simulation experiment data #
# load("arco.Rdata")
# load("arcodrift1.RData")
# load("arcodrift2.RData")
# load("diff_arco.RData")
# load("diff_arcodrift1.RData")
# load("diff_arcodrift2.RData")
# load("sc.Rdata")
# load("scdrift1.RData")
# load("scdrift2.RData")
# load("diff_sc.RData")
# load("diff_scdrift1.RData")
# load("diff_scdrift2.RData")

##############################################################################
############################## 1. set variables ##############################  
##############################################################################

# input for the data generating function: random walk process
y0 = 0 # initial value
delta1 = -0.5 # negative drift
delta2 = 0.5 # positive drift
variance = 1 # variance of innovations

# counterfactual experiment information
N = 20 # amount of units
T = 100 # length time series for each unit
t0 = T/2 # time of "intervention"
tprior <- t0-1 # time before "intervention"
R = 500 # Number of bootstrap replications in the ArCo to construct confidence intervals
Y = 1000 # experiment size, per data category

##############################################################################
############### 3. run function & generate experiment results ################
##############################################################################

# model = "ArCo" #

ptm <- proc.time() # start ArCo timer

# 'normal' RW data - diff = FALSE #

  # no drift #
  arco<- replicate(n=Y, nonstationary_ts_counterfactual("ArCo", T, N, y0, 0, variance, R, FALSE))
  save(arco, file="arco.RData")
  # negative drift #
  arcodrift1<- replicate(n=Y, nonstationary_ts_counterfactual("ArCo", T, N, y0, delta1, variance, R, FALSE))
  save(arcodrift1, file="arcodrift1.RData")
  # positive drift #
  arcodrift2<- replicate(n=Y, nonstationary_ts_counterfactual("ArCo", T, N, y0, delta2, variance, R, FALSE))
  save(arcodrift2, file="arcodrift2.RData")

# differenced RW data - diff = TRUE #

  # no drift #
  diff_arco<- replicate(n=Y, nonstationary_ts_counterfactual("ArCo", T, N, y0, 0, variance, R, TRUE))
  save(diff_arco, file="diff_arco.RData")
  # negative drift #
  diff_arcodrift1<- replicate(n=Y, nonstationary_ts_counterfactual("ArCo", T, N, y0, delta1, variance, R, TRUE))
  save(diff_arcodrift1, file="diff_arcodrift1.RData")
  # positive drift #
  diff_arcodrift2<- replicate(n=Y, nonstationary_ts_counterfactual("ArCo", T, N, y0, delta2, variance, R, TRUE))
  save(diff_arcodrift2, file="diff_arcodrift2.RData")

time_arco <- proc.time() - ptm # stop ArCo timer

# model = "SC" #

ptm <- proc.time() # start SC timer

# 'normal' RW data - diff = FALSE #

  # no drift #
  sc <- replicate(n=Y, nonstationary_ts_counterfactual("SC", T, N, y0, 0, variance, R, FALSE))
  save(sc, file="sc.RData")
  # negative drift #
  scdrift1 <- replicate(n=Y, nonstationary_ts_counterfactual("SC", T, N, y0, delta1, variance, R, FALSE))
  save(scdrift1, file="scdrift1.RData")
  # positive drift #
  scdrift2 <- replicate(n=Y, nonstationary_ts_counterfactual("SC", T, N, y0, delta2, variance, R, FALSE))
  save(scdrift2, file="scdrift2.RData")

# differenced RW data - diff = TRUE #

  # no drift #
  diff_sc <- replicate(n=Y, nonstationary_ts_counterfactual("SC", T, N, y0, 0, variance, R, TRUE))
  save(diff_sc, file="diff_sc.RData")
  # negative drift #
  diff_scdrift1 <- replicate(n=Y, nonstationary_ts_counterfactual("SC", T, N, y0, delta1, variance, R, TRUE))
  save(diff_scdrift1, file="diff_scdrift1.RData")
  # positive drift #
  diff_scdrift2 <- replicate(n=Y, nonstationary_ts_counterfactual("SC", T, N, y0, delta2, variance, R, TRUE))
  save(diff_scdrift2, file="drift_scdrift2.RData")

time_sc <- proc.time() - ptm # stop SC timer

##############################################################################
################################# 4. results #################################  
##############################################################################

# percentage of significant ArCo results on (differenced) Random Walks, for Y=1000 per category # 

alpha <- c(0.001, 0.01, 0.05, 0.1)

test0 <- percentagesignificant(arco, alpha)
test1 <- percentagesignificant(arcodrift1, alpha)
test2 <- percentagesignificant(arcodrift2, alpha)
tests <- rbind(test0,test1,test2)
colnames(tests)<- alpha

diff_test0 <- percentagesignificant(diff_arco, alpha)
diff_test1 <- percentagesignificant(diff_arcodrift1, alpha)
diff_test2 <- percentagesignificant(diff_arcodrift2, alpha)
diff_tests <- rbind(diff_test0,diff_test1,diff_test2)
colnames(diff_tests)<- alpha

rm(test0,test1,test2,diff_test0,diff_test1,diff_test2)

rbind(tests, diff_tests) # probability of rejections for normal data and first difference data
colSums(tests)/3 # average probability of rejections for normal data
colSums(diff_tests)/3 # average probability of rejections for first diff data

# extract ArCo delta's #
delta_arco<- arco["delta",]
delta_arcodrift1<- arcodrift1["delta",]
delta_arcodrift2<- arcodrift2["delta",]
delta_diff_arco<- diff_arco["delta",]
delta_diff_arcodrift1<- diff_arcodrift1["delta",]
delta_diff_arcodrift2<- diff_arcodrift2["delta",]

# compare ArCo & SC magnitudes by comparing gaps (gaps w.r.t. original simulated data Y1) #

# construct gaps per arco
gaps_arco <- arcogaps(arco, diff=FALSE)
gaps_arcodrift1 <- arcogaps(arcodrift1, diff=FALSE)
gaps_arcodrift2 <- arcogaps(arcodrift2, diff=FALSE)
gaps_diff_arco <- arcogaps(diff_arco, diff=TRUE)
gaps_diff_arcodrift1 <- arcogaps(diff_arcodrift1, diff=TRUE)
gaps_diff_arcodrift2 <- arcogaps(diff_arcodrift2, diff=TRUE)

# sc
# construct gaps per sc
gaps_sc <- as.data.frame(sc["gaps",])
gaps_scdrift1 <- as.data.frame(scdrift1["gaps",])
gaps_scdrift2 <- as.data.frame(scdrift2["gaps",])
gaps_diff_sc <- as.data.frame(diff_sc["gaps",])
gaps_diff_scdrift1 <- as.data.frame(diff_scdrift1["gaps",])
gaps_diff_scdrift2 <- as.data.frame(diff_scdrift2["gaps",])

# take sums over all gaps
# for each category within arco and SC: sum over T, for all Y=1000 simulations
# arco 
sum_gaps_arco <- colSums(gaps_arco) # each gap is small delta
sum_gaps_arcodrift1 <- colSums(gaps_arcodrift1)
sum_gaps_arcodrift2 <- colSums(gaps_arcodrift2)
sum_gaps_diff_arco <- colSums(gaps_diff_arco) 
sum_gaps_diff_arcodrift1 <- colSums(gaps_diff_arcodrift1)
sum_gaps_diff_arcodrift2 <- colSums(gaps_diff_arcodrift2)

# sc
sum_gaps_sc <- colSums(gaps_sc) 
sum_gaps_scdrift1 <- colSums(gaps_scdrift1)
sum_gaps_scdrift2 <- colSums(gaps_scdrift2)
sum_gaps_diff_sc <- colSums(gaps_diff_sc) 
sum_gaps_diff_scdrift1 <- colSums(gaps_diff_scdrift1)
sum_gaps_diff_scdrift2 <- colSums(gaps_diff_scdrift2)

# plot sums of gaps (= total surface of gaps)
par(mfrow = c(3, 2))  # 3 rows and 2 columns  

# arco
plot(sum_gaps_arco)
plot(sum_gaps_diff_arco)
plot(sum_gaps_arcodrift1)
plot(sum_gaps_diff_arcodrift1)
plot(sum_gaps_arcodrift2)
plot(sum_gaps_diff_arcodrift2)

hist(sum_gaps_arco)
hist(sum_gaps_diff_arco)
hist(sum_gaps_arcodrift1)
hist(sum_gaps_diff_arcodrift1)
hist(sum_gaps_arcodrift2)
hist(sum_gaps_diff_arcodrift2)

# sc
plot(sum_gaps_sc)
plot(sum_gaps_diff_sc)
plot(sum_gaps_scdrift1)
plot(sum_gaps_diff_scdrift1)
plot(sum_gaps_scdrift2)
plot(sum_gaps_diff_scdrift2)

hist(sum_gaps_sc)
hist(sum_gaps_diff_sc)
hist(sum_gaps_scdrift1)
hist(sum_gaps_diff_scdrift1)
hist(sum_gaps_scdrift2)
hist(sum_gaps_diff_scdrift2)

# sum difference per category
# the mean reported in the t-tests below is NOT the same as the mean ATET (=Delta_t) over all Y's 
# this will be performed below for completeness
# however: the results (t-tests and p-values) are the same, the means are just different.
# the mean reported in the t-tests below is the mean over all Y's of the SUM of all gaps  (=sum(delta_t)) 

# 1. simple RW
t1 <- t.test(sum_gaps_arco, sum_gaps_sc)
# 2. RW + negative drift
t2 <- t.test(sum_gaps_arcodrift1, sum_gaps_scdrift1)
# 3. RW + positive drift 
t3 <- t.test(sum_gaps_arcodrift2, sum_gaps_scdrift2)
# 4. simple RW, first difference
t4 <- t.test(sum_gaps_diff_arco, sum_gaps_diff_sc)
# 5. RW + negative drift, first difference
t5 <- t.test(sum_gaps_diff_arcodrift1, sum_gaps_diff_scdrift1)
# 6. RW + positive drift, first difference    
t6 <- t.test(sum_gaps_diff_arcodrift2, sum_gaps_diff_scdrift2) 

#bind tests in a table
t.tests.sum<-as.data.frame(rbind(t1, t2, t3, t4, t5, t6))
# write_xlsx(t.tests.sum, 't.tests.SUM.xlsx')


# ATET difference per category # 
# ATET=Delta_t= 1/(T-T0+1) sum_{t=T0}^T delta_t, where delta_t is the gap at each t.    
# the mean reported in the t-tests below is the mean over all Y's all Delta_t's

# ATET per data category & per simulation
# arco 
Delta_t_arco <- sum_gaps_arco/(T-t0+1) # ATET per simulation # note how this is the same as function delta
Delta_t_arcodrift1<-sum_gaps_arcodrift1/(T-t0+1)
Delta_t_arcodrift2<-sum_gaps_arcodrift2/(T-t0+1)

Delta_t_diff_arco <- sum_gaps_diff_arco/(T-t0) # -1 because of of first difference 
Delta_t_diff_arcodrift1<-sum_gaps_diff_arcodrift1/(T-t0)
Delta_t_diff_arcodrift2<-sum_gaps_diff_arcodrift2/(T-t0)

# sc
Delta_t_sc <- sum_gaps_sc/(T-t0+1) # ATET per simulation
Delta_t_scdrift1<-sum_gaps_scdrift1/(T-t0+1)
Delta_t_scdrift2<-sum_gaps_scdrift2/(T-t0+1)

Delta_t_diff_sc <- sum_gaps_diff_sc/(T-t0) # -1 because of of first difference 
Delta_t_diff_scdrift1<-sum_gaps_diff_scdrift1/(T-t0)
Delta_t_diff_scdrift2<-sum_gaps_diff_scdrift2/(T-t0)

# 1. simple RW
t1.ATET <- t.test(Delta_t_arco, Delta_t_sc)
# 2. RW + negative drift
t2.ATET <- t.test(Delta_t_arcodrift1, Delta_t_scdrift1)
# 3. RW + positive drift 
t3.ATET <- t.test(Delta_t_arcodrift2, Delta_t_scdrift2)
# 4. simple RW, first difference
t4.ATET <- t.test(Delta_t_diff_arco, Delta_t_diff_sc)
# 5. RW + negative drift, first difference
t5.ATET <- t.test(Delta_t_diff_arcodrift1, Delta_t_diff_scdrift1)
# 6. RW + positive drift, first difference    
t6.ATET <- t.test(Delta_t_diff_arcodrift2, Delta_t_diff_scdrift2)

# bind tests in a table
t.tests.ATET<-as.data.frame(rbind(t1.ATET, t2.ATET, t3.ATET, t4.ATET, t5.ATET, t6.ATET))
# write_xlsx(t.tests.ATET, 't.tests.ATET.xlsx')
