### ### ### ### ### SIMULATION EXPERIMENT ### ### ### ### ### 

# 14 juni 2021 
# v4 
# M Hagens

setwd("/cloud/project")
cat("\014")        
rm(list = ls())

source("randomwalk.R")
library(ArCo)
library(Synth)
library(glmnet) # contains LASSO
library(reshape2) # contains melt(), when using SC 

set.seed(42)

#### set variables #### 
# alpha=0.05
y0=0
delta=0 
#variance=(0.5^2)
variance = 0.5

#n=29
#t=29
N=29
T=100 

t0=T/2
tprior <- t0-1
R = 1000 # Number of bootstrap replications in the arco
# diff = FALSE
# model = ArCo #model can either be "ArCo" or "SC"

counterfactual_NS_exp <- function(model, T, N, y0, delta, variance, R, diff){
  
  #### using random walk function #### 
  RW <- ts(replicate(n = N, # number of units
                     randomwalk(T, y0, delta, variance))) 
  
  ####  if  specified, difference the data   #### 
  if(diff == TRUE){
    RW <-diff(RW) 
    T <- nrow(RW)
  }
  
  #### experiment  ####
  randomtreated <- sample(1:N, 1)
  
  if(model == "ArCo"){ # put dataframe in ArCo format, run and return ArCo
  
  list<-list(RW)
  arco <- fitArCo(data = list,  fn = cv.glmnet, p.fn = predict, 
                  treated.unit = randomtreated, t0 = t0, 
                  boot.cf=TRUE, R = R, # should I change the bootstrap?
                  VCOV.type = "iid"
  )
  
  # 1. open jpeg file
  #jpeg("plots/rplot.jpg")
  # 2. create the plot
  # plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Random Walk ArCo")
  # 3. close the file
  #dev.off()
  
  return(arco)
  ls()
  
  }
  
  if(model == "SC"){ # put dataframe in SC format, run and return SC

    # the melt() function can reshape a long data frame to wide data frame, 
    # which is needed as input for SC
    dfLong <- melt(RW, value.name="X")
    colnames(dfLong) <- c("Time","Unit","X")
    
    # create a numeric unit variable to indicate which unit is which "unit.num"
    n <- ncol(RW)
    t <- nrow(RW)
    dfLong$unit.num <- rep(c(1:n),each=t)
    #dfLong$unit.num <- rep(c(1:N), each = T)
    
    dfLong$Unit <- as.character(dfLong$Unit)
    
    donorpool <- rep(c(1:N))
    donorpool <- donorpool[c(-randomtreated)]
    
    # create matrix from panel data that provides input for synth()
    dataprep.out<-
      dataprep(
        foo = dfLong,
        predictors = "X",
        predictors.op = "mean", # mean is the default
        dependent = "X",
        unit.variable = "unit.num",  # this has to be a numeric variable
        time.variable = "Time",
        treatment.identifier = randomtreated,
        controls.identifier = donorpool,
        time.predictors.prior = c(1:tprior),
        time.optimize.ssr = c(1:t0), 
        unit.names.variable = "Unit", 
        time.plot = 1:T
      )
    
    # run synth
    synth <- synth(dataprep.out, method = "BFGS") #why BFGS?
    #X0 the control cases after the treatment
    #X1 the control case before the treatment
    #Z1 the treatment case before the treatment
    #Z0 the treatment case after the treatment
    
    # calculate difference in trend between simulated data and its synthetic control
    gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth$solution.w)
    #print(gaps[1:3, 1]) # check
    
    # pre built tables from synth objects
    synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth)
    names(synth.tables) # check
    
    # comparing pre-treatment predictor values for the treated unit, the synthetic control unit, and all the units in the sample
    
    # synth.tables$tab.pred # check balance across treated and control for pre-period predictors
    
    # plot treatment vs control outcomes for pre and post periods
    path.plot(synth.res = synth, dataprep.res = dataprep.out,
              Ylab = "Value", Xlab = "Time",
              Legend = c("Random Walk","Synthetic Control"), 
              Legend.position = "bottomright",
              Main = "Simulated RW data and its SC")
    
    # gap plot
    gaps.plot(synth.res = synth, dataprep.res = dataprep.out,
              Ylab = "Value", Xlab = "Time",
              Main = "Gap between RW and SC")
    
    return(synth)
    ls()
  }
}

ptm <- proc.time() # start timer

# run function & generate experiment results
#ArCo
Xarco<- replicate(n=2, counterfactual_NS_exp("ArCo", T, N, y0, delta, variance, R, TRUE))
Yarco<- replicate(n=2, counterfactual_NS_exp("ArCo", T, N, y0, delta, variance, R, FALSE))
#SC
Xsc<- replicate(n=2, counterfactual_NS_exp("SC", T, N, y0, delta, variance, R, TRUE))
Ysc<- replicate(n=2, counterfactual_NS_exp("SC", T, N, y0, delta, variance, R,FALSE))

time <- proc.time() - ptm # stop timer

#### significance tests ####
alpha <- c(0.001, 0.01, 0.05, 0.1)
Xarco_test_0.001 <- (Xarco["p.value",]   < alpha[1])
Xarco_test_0.01 <- (Xarco["p.value",]   < alpha[2])
Xarco_test_0.05 <- (Xarco["p.value",]   < alpha[3])
Xarco_test_0.10 <- (Xarco["p.value",]   < alpha[4])
Xarco_test<- cbind(Xarco_test_0.001, Xarco_test_0.01, Xarco_test_0.05,  Xarco_test_0.10)
Yarco_test_0.001 <- (Yarco["p.value",]   < alpha[1])
Yarco_test_0.01 <- (Yarco["p.value",]   < alpha[2])
Yarco_test_0.05 <- (Yarco["p.value",]   < alpha[3])
Yarco_test_0.10 <- (Yarco["p.value",]   < alpha[4])
Yarco_test<- cbind(Yarco_test_0.001, Yarco_test_0.01, Yarco_test_0.05,  Yarco_test_0.10)


### SC does not have any p values
# Xsc_test_0.001 <- (Xsc["p.value",]   < alpha[1])
# Xsc_test_0.01 <- (Xsc["p.value",]   < alpha[2])
# Xsc_test_0.05 <- (Xsc["p.value",]   < alpha[3])
# Xsc_test_0.10 <- (Xsc["p.value",]   < alpha[4])
# Xsc_test<- cbind(Xsc_test_0.001, Xsc_test_0.01, Xsc_test_0.05,  Xsc_test_0.10)
# Ysc_test_0.001 <- (Ysc["p.value",]   < alpha[1])
# Ysc_test_0.01 <- (Ysc["p.value",]   < alpha[2])
# Ysc_test_0.05 <- (Ysc["p.value",]   < alpha[3])
# Ysc_test_0.10 <- (Ysc["p.value",]   < alpha[4])
# Ysc_test<- cbind(Ysc_test_0.001, Ysc_test_0.01, Ysc_test_0.05,  Ysc_test_0.10)
