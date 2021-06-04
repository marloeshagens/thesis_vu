### ### ### ### ### SIMULATION EXPERIMENT ### ### ### ### ### 

# 4 juni 2021 
# v3 
# M Hagens

setwd("/cloud/project")
cat("\014")        
rm(list = ls())

source("randomwalk.R")
set.seed(123)
library(ArCo) # contains ArCo package
library(glmnet) # contains LASSO

# simulate per n = 29, T = 29 ...
n=29
t=29
N=100
T=1000

# initialization innovation u_t
var = 1 #variance of normally distributed error
drift = -0.5 # this can be a drift. if 0, then it is unbiased; equally likely to go up or down.  

#### using ARIMA (0,1,0) #### 
# RW1. a random walk (simple and independently)
RW1 <- ts(replicate(n = 5, # number of units (paths to generate) 
                    arima.sim(model = list(order = c(0, 1, 0)),
                              n = t, # number of observations
                              mean = 0,
                              sd = sqrt(var)
                    )))
matplot(RW1, 
        type = "l", 
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"), 
        lty = 1, 
        lwd = 2,
        main = "Random Walk (1), no drift: mean = 0, sd=1",
        xlab = "Time",
        ylab = "Value")

# RW2. a random walk + drift
RW2 <- ts(replicate(n = 5, # number of units (paths to generate) 
                    arima.sim(model = list(order = c(0, 1, 0)),
                              n = t, # number of observations
                              mean = -0.5,
                              sd = 1,
                    )))
matplot(RW2, 
        type = "l", 
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"), 
        lty = 1, 
        lwd = 2,
        main = "Random walk (2), with drift: mean =-0.5, sd=1",
        xlab = "Time",
        ylab = "Value")

# RW3. a random walk + trend
trend <- -0.5
RW3 <- ts(replicate(n = 5, # number of units (paths to generate)
                    trend + arima.sim(model = list(order = c(0, 1, 0)),
                              n = t, # number of observations
                              mean = -0.5,
                              sd = 1,
                    )))

matplot(RW3,
        type = "l",
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"),
        lty = 1,
        lwd = 2,
        main = "Random walk (3), with drift + trend: mean =-0.5, sd=1",
        xlab = "Time",
        ylab = "Value")

#### using randomwalk function #### 
#RW1_ <- randomwalk(T, 0, 0, 1) 
RW1_ <- ts(replicate(n = N, # number of units
                     randomwalk(T, 0, 0, 1))) # without drift
RW2_ <- ts(replicate(n = N, 
                     randomwalk (T, 0, drift, 1))) # with drift - -0.5
matplot(RW1_,
        type = "l",
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"),
        lty = 1,
        lwd = 2,
        main = "Random walk",
        xlab = "Time",
        ylab = "Value")
matplot(RW2_,
        type = "l",
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"),
        lty = 1,
        lwd = 2,
        main = "Random walk + drift",
        xlab = "Time",
        ylab = "Value")

# differenced data

# DiffRW1
RW_diff <-diff(RW1_)
matplot(RW_diff,
        type = "l",
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"),
        lty = 1,
        lwd = 2,
        main = "First Order Difference (RW)",
        xlab = "Time",
        ylab = "Diff Value")

# DiffRW2
RW_drift_diff <-diff(RW2_)
matplot(RW_drift_diff,
        type = "l",
        col = c("steelblue", "darkgreen", "darkred", "orange", "pink"),
        lty = 1,
        lwd = 2,
        main = "First Order Difference (RW with Drift)",
        xlab = "Time",
        ylab = "Diff Value")

#### experiment example #### 
randomtreated <- sample(1:N,1)
t0 <- T/2
list<-list(RW1_)# dataprep arco
arco <- fitArCo(data = list,  fn = cv.glmnet, p.fn = predict, 
                 treated.unit = randomtreated, t0 = t0, 
                 boot.cf=TRUE, R=1000, # should I change the bootstrap?
                 VCOV.type = "iid"
)
plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Standard Brownian Motion ArCo")
arco$p.value  < 0.05

