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
diff = FALSE

counterfactual_NS_arco_experiment <- function(T, y0, delta, variance, R, diff){
  
  #### using randomwalk function #### 
  RW <- ts(replicate(n = n, # number of units
                     randomwalk(T, y0, delta, variance))) 
  
  if(diff == TRUE){
    RW <-diff(RW) # difference the data
  }
  
  list<-list(RW) # dataprep arco
  
  #### experiment  #### 
  randomtreated <- sample(1:n, 1)
  arco <- fitArCo(data = list,  fn = cv.glmnet, p.fn = predict, 
                  treated.unit = randomtreated, t0 = t0, 
                  boot.cf=TRUE, R=1000, # should I change the bootstrap?
                  VCOV.type = "iid"
  )
  
  # 1. open jpeg file
  #jpeg("plots/rplot.jpg")
  # 2. create the plot
  plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Random Walk ArCo")
  # 3. close the file
  #dev.off()
  return(arco)
  ls()
}

# Start the clock!
ptm <- proc.time()

# run function
X<- replicate(n=2, counterfactual_NS_arco_experiment(T, y0, delta, variance,R,TRUE))

# Stop the clock
proc.time() - ptm

Y<- replicate(n=2, counterfactual_NS_arco_experiment(T, y0, delta, variance,R,FALSE))

# create dataframe of test values
alpha <- c(0.001, 0.01, 0.05, 0.1)
Xtest_0.001 <- (X["p.value",]   < alpha[1])
Xtest_0.01 <- (X["p.value",]   < alpha[2])
Xtest_0.05 <- (X["p.value",]   < alpha[3])
Xtest_0.10 <- (X["p.value",]   < alpha[4])
Xtest<- cbind(Xtest_0.001, Xtest_0.01, Xtest_0.05,  Xtest_0.10 )
Ytest_0.001 <- (X["p.value",]   < alpha[1])
Ytest_0.01 <- (X["p.value",]   < alpha[2])
Ytest_0.05 <- (X["p.value",]   < alpha[3])
Ytest_0.10 <- (X["p.value",]   < alpha[4])
Ytest<- cbind(Ytest_0.001, Ytest_0.01, Ytest_0.05,  Ytest_0.10 )