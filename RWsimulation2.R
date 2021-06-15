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
t0 <- T/2

# amount of simulations
sim = 2

# initialization innovation u_t
var = 1 #variance of normally distributed error
drift = -0.5 # this can be a drift. if 0, then it is unbiased; equally likely to go up or down.  

# maak hier een functie van ; return test
test<- rep(NA, sim)
set.seed(123)
for (i in 1:sim) {
  #### using randomwalk function #### 
  RW1 <- ts(replicate(n = n, # number of units
                      randomwalk(T, 0, 0, 1))) # without drift
  RW2 <- ts(replicate(n = n, 
                      randomwalk (T, 0, drift, 1))) # with drift - -0.5
  
  # differenced data
  RW1_diff <-diff(RW1)
  RW2_diff <-diff(RW2)
  
  #### experiment  #### 
  randomtreated <- sample(1:n,1)
  list<-list(RW1)# dataprep arco
  arco <- fitArCo(data = list,  fn = cv.glmnet, p.fn = predict, 
                  treated.unit = randomtreated, t0 = t0, 
                  boot.cf=TRUE, R=10, # should I change the bootstrap?
                  VCOV.type = "iid"
  )
  #plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Standard Brownian Motion ArCo")
  test[i] <- arco$p.value  < 0.05
}


counterfactual_NS_experiment <- function(sim, alpha){
  set.seed(123)
  test1<- rep(NA, sim)
  test2<- rep(NA, sim)
  #delta<- rep(NA, sim)
  for (i in 1:sim) {
    #### using randomwalk function #### 
    RW1 <- ts(replicate(n = n, # number of units
                        randomwalk(T, 0, 0, 1))) # without drift
    RW2 <- ts(replicate(n = n, 
                        randomwalk (T, 0, drift, 1))) # with drift - -0.5
    
    # differenced data
    RW1_diff <-diff(RW1)
    RW2_diff <-diff(RW2)
    
    #### experiment  #### 
    randomtreated <- sample(1:n,1)
    list1<-list(RW1)# dataprep arco
    arco1 <- fitArCo(data = list1,  fn = cv.glmnet, p.fn = predict, 
                    treated.unit = randomtreated, t0 = t0, 
                    boot.cf=TRUE, R=100, # should I change the bootstrap?
                    VCOV.type = "iid"
    )
    #plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Standard Brownian Motion ArCo")
    test1[i] <- arco1$p.value  < alpha
    
    list2<-list(RW2)# dataprep arco
    arco2 <- fitArCo(data = list2,  fn = cv.glmnet, p.fn = predict, 
                     treated.unit = randomtreated, t0 = t0, 
                     boot.cf=TRUE, R=100, # should I change the bootstrap?
                     VCOV.type = "iid"
    )
    #plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Standard Brownian Motion ArCo")
    test2[i] <- arco1$p.value  < alpha
    }
  return(test1)
  return(test2)
}

sim=5
alpha=0.05
test<- counterfactual_NS_experiment(sim,alpha)