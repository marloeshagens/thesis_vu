### ### ### ### ### SIMULATION EXPERIMENT ### ### ### ### ### 

# 14 juni 2021 
# v4 
# M Hagens

setwd("/cloud/project")
cat("\014")        
rm(list = ls())

source("randomwalk.R")
library(ArCo) # contains ArCo package
library(glmnet) # contains LASSO

set.seed(42)

ptm <- proc.time() # Start the clock!

#### set variables #### 
alpha=0.05
y0=0
delta=0 
variance=(0.5^2)
source("randomwalk.R")

n=29
t=29
N=100
T=100 

t0=T/2
R = 1000 # Number of bootstrap replications in the arco

counterfactual_NS_arco_experiment <- function(T, y0, delta, variance, R){
  
  #### using randomwalk function #### 
  RW <- ts(replicate(n = n, # number of units
                        randomwalk(T, y0, delta, variance))) 

    # differenced data
    #RW1_diff <-diff(RW1)
    #RW2_diff <-diff(RW2)
    
    #### experiment  #### 
    randomtreated <- sample(1:n, 1)
    
    list<-list(RW) # dataprep arco
    
    arco <- fitArCo(data = list,  fn = cv.glmnet, p.fn = predict, 
                     treated.unit = randomtreated, t0 = t0, 
                     boot.cf=TRUE, R=1000, # should I change the bootstrap?
                     VCOV.type = "iid"
    )
    plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Random Walk ArCo")
    return(arco)
    ls()
}
X<- replicate(n=5, counterfactual_NS_arco_experiment(T, y0, delta, variance,R))

alpha <- c(0.001, 0.01, 0.05, 0.1)
test_0.001 <- (X["p.value",]   < alpha[1])
test_0.01 <- (X["p.value",]   < alpha[2])
test_0.05 <- (X["p.value",]   < alpha[3])
test_0.10 <- (X["p.value",]   < alpha[4])
test<- cbind(test_0.001, test_0.01, test_0.05,  test_0.10 )

proc.time() - ptm # Stop the clock

