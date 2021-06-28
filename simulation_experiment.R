### ### ### ### ### SIMULATION EXPERIMENT ### ### ### ### ### 

# 24 juni 2021 
# v6
# author: M Hagens

setwd("~/Documents/thesis_vu-master")
cat("\014")        
rm(list = ls())

##############################################
############# 0. load functions ##############
##############################################

source("randomwalk.R") # data generating process function #
library(ArCo) # ArCo package #
library(Synth) # Synthetic Control package #
library(glmnet) # contains LASSO #
library(reshape2) # contains melt(), when using SC #
library(Matrix) # used for colSums #

set.seed(42)

##############################################
############# 1. set variables ############### 
##############################################
y0=0
delta1=-0.5
delta2=0.5
variance = 1

N=20
T=100 

t0=T/2
tprior <- t0-1
R = 500 # Number of bootstrap replications in the arco

##############################################
# 2. define count.fact. experiment function  #
##############################################

nonstationary_experiment <- function(model, T, N, y0, delta, variance, R, diff){
# this function creates a N*T matrix of random walk data, and puts in into a counterfactual method model # 
# model: either "ArCo" or "SC" #
# T:  the length of each time series #
# N: the number of units #
# delta: the (possible) drift in the random walk # 
# variance: specifies the random walk's error term distribution #
# R: if model = "ArCo", R specifies the number of bootstrap replications in the ArCo #
  
  # DGP using random walk function # 
  RW <- ts(replicate(n = N, # number of units 
                     randomwalk(T, y0, delta, variance))) # produces random walk per unit
  
  # if specified, difference the data # 
  if(diff == TRUE){
    RW <-diff(RW) 
    T <- nrow(RW)
  }
  
  # experiment #
  randomtreated <- sample(1:N, 1)
  
    if(model == "ArCo"){ # put dataframe in ArCo format, run and return ArCo
      list<-list(RW)
      
      arco <- fitArCo(data = list,  fn = cv.glmnet, p.fn = predict, 
                      treated.unit = randomtreated, t0 = t0, 
                      boot.cf=TRUE, R = R,
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
  
    if(model == "SC"){ # put dataframe in SC format, run and return SC # 

      dfLong <- melt(RW, value.name="X") # reshapes data frame for the SC input # 
      colnames(dfLong) <- c("Time","Unit","X")
      
      n <- ncol(RW)
      t <- nrow(RW)
      dfLong$unit.num <- rep(c(1:n),each=t) # create a numeric unit variable to indicate which unit is which # 
      
      dfLong$Unit <- as.character(dfLong$Unit) # Unit column should be made of characters in order for the SC to work # 
      
      donorpool <- rep(c(1:N)) # define donorpool # 
      donorpool <- donorpool[c(-randomtreated)]
      
      # create matrix from panel data that provides input for synth()
      dataprep.out<-
        dataprep(
          foo = dfLong,
          predictors = "X",
          predictors.op = "mean", # mean is the default
          dependent = "X",
          unit.variable = "unit.num",
          time.variable = "Time",
          treatment.identifier = randomtreated,
          controls.identifier = donorpool,
          time.predictors.prior = c(1:tprior),
          time.optimize.ssr = c(1:t0), 
          unit.names.variable = "Unit", 
          time.plot = 1:T
        )
      
      # run synth
      synth <- synth(dataprep.out, method = "BFGS") # why BFGS?
      #X0 the control cases after the treatment
      #X1 the control case before the treatment
      #Z1 the treatment case before the treatment
      #Z0 the treatment case after the treatment
      
      # calculate difference in trend between simulated data and its synthetic control
      gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth$solution.w)
      #print(gaps[1:3, 1]) 
      
      # pre built tables from synth objects
      synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth)
      names(synth.tables) 

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

##############################################
# 3. run func. & generate experiment results #
##############################################

# experiment sample size, per category #
Y=1000

############### model = "ArCo" ############### 

ptm <- proc.time() # start timer

# 'normal' RW data - diff = FALSE #

  # no drift #
  arco<- replicate(n=1000, nonstationary_experiment("ArCo", T, N, y0, 0, variance, R, FALSE))
  time_arco <- proc.time() - ptm # stop timer
  # negative drift #
  arcodrift1<- replicate(n=Y, nonstationary_experiment("ArCo", T, N, y0, delta1, variance, R, FALSE))
  # positive drift #
  arcodrift2<- replicate(n=Y, nonstationary_experiment("ArCo", T, N, y0, delta2, variance, R, FALSE))

# differenced RW data - diff = TRUE #

  # no drift #
  diff_arco<- replicate(n=Y, nonstationary_experiment("ArCo", T, N, y0, 0, variance, R, TRUE))
  # negative drift #
  diff_arcodrift1<- replicate(n=Y, nonstationary_experiment("ArCo", T, N, y0, delta1, variance, R, TRUE))
  # positive drift #
  diff_arcodrift2<- replicate(n=Y, nonstationary_experiment("ArCo", T, N, y0, delta2, variance, R, TRUE))

# percentage of significant ArCo results on (differenced) Random Walks, for Y=1000 per category # 
  
alpha <- c(0.001, 0.01, 0.05, 0.1)

percentagesignificant <- function(arcooutput, alpha){
  test <- matrix(NA, length(arcooutput["p.value",]), length(alpha))
  for (i in 1:length(alpha)) {
    test[,i] <- (arcooutput["p.value",] < alpha[i])
  }
  percentage_true <- colSums(test == TRUE)/nrow(test)
  return(percentage_true)
}

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

rbind(tests, diff_tests)

time_arco <- proc.time() - ptm # stop timer

################ model = "SC" ################ 

ptm <- proc.time()

# 'normal' RW data - diff = FALSE #

  # no drift #
  sc<- replicate(n=Y, nonstationary_experiment("SC", T, N, y0, 0, variance, R, FALSE))
  # negative drift #
  scdrift1<- replicate(n=Y, nonstationary_experiment("SC", T, N, y0, delta1, variance, R, FALSE))
  # positive drift #
  scdrift2<- replicate(n=Y, nonstationary_experiment("SC", T, N, y0, delta2, variance, R, FALSE))

# differenced RW data - diff = TRUE #

  # no drift #
  drift_sc<- replicate(n=Y, nonstationary_experiment("SC", T, N, y0, 0, variance, R, TRUE))
  # negative drift #
  drift_scdrift1<- replicate(n=Y, nonstationary_experiment("SC", T, N, y0, delta1, variance, R, TRUE))
  # positive drift #
  drift_scdrift2<- replicate(n=Y, nonstationary_experiment("SC", T, N, y0, delta2, variance, R, TRUE))

time_sc <- proc.time() - ptm # stop timer

##############################################
###### 4. compare ArCo & SC magnitudes #######
##############################################

#extract delta's
arco_delta<- arco["delta",]
arcodrift1_delta<- arcodrift1["delta",]
arcodrift2_delta<- arcodrift2["delta",]

diff_arco_delta<- diff_arco["delta",]
diff_arcodrift1_delta<- diff_arcodrift1["delta",]
diff_arcodrift2_delta<- diff_arcodrift2["delta",]

# drift_sc_delta<- drift_sc["delta",]
# drift_scdrift1_delta<- drift_scdrift1["delta",]
# drift_scdrift2_delta<- drift_scdrift2["delta",]
# sc_delta<- sc["delta",]
# scdrift1_delta<- scdrift1["delta",]
# scdrift2_delta<- scdrift2["delta",]


# pick a random dataset and plot the ArCo and the SC #
# maybe put this in an entirely different code? at least, at the bottom
# works for the non-differenced way 
synthcf=dataprep.out$Y0plot%*%synth$solution.w
plot(dataprep.out$Y1plot,type="l",ylab="Y",xlab="Time")
lines(synthcf,col=2)
abline(v=50,col=4,lty=2)
lines(c(fitted(arco)[,1], arco$cf[,1]),col=4)
legend("bottomleft",legend=c("Y","Synth","ArCo"),
       col=c(1,2,4),
       cex = 1,
       seg.len = 1,bty = "n",lty=1)

# save results 
# write.csv(arco,'arco.csv')
# write.csv(arcodrift1,'arcodrift1.csv')
# write.csv(arcodrift2,'arcodrift2.csv')


# prop.test(80, 100, p = 0, alternative = "two.sided", correct = TRUE)
# prop.test(85, 100, p = 0, alternative = "two.sided", correct = TRUE)
# prop.test(88, 100, p = 0, alternative = "two.sided", correct = TRUE)
# prop.test(89, 100, p = 0, alternative = "two.sided", correct = TRUE)
# 
# prop.test(79, 100, p = 0, alternative = "two.sided", correct = TRUE)
# prop.test(84, 100, p = 0, alternative = "two.sided", correct = TRUE)
# prop.test(90, 100, p = 0, alternative = "two.sided", correct = TRUE)
# prop.test(90, 100, p = 0, alternative = "two.sided", correct = TRUE)
# 
# prop.test(0, 100, p = 0.1, alternative = "two.sided", correct = TRUE)
# prop.test(2, 100, p = 0.1, alternative = "two.sided", correct = TRUE)
# prop.test(4, 100, p = 0.1, alternative = "two.sided", correct = TRUE)
# prop.test(9, 100, p = 0.1, alternative = "two.sided", correct = TRUE)


# #### old notation significance tests ####
# arco_test_0.001 <- (arco["p.value",]   < alpha[1])
# arco_test_0.01 <- (arco["p.value",]   < alpha[2])
# arco_test_0.05 <- (arco["p.value",]   < alpha[3])
# arco_test_0.10 <- (arco["p.value",]   < alpha[4])
# arco_test<- cbind(arco_test_0.001, arco_test_0.01, arco_test_0.05,  arco_test_0.10)
# test0<- colSums(arco_test == TRUE)/nrow(arco_test)*100 #colSums(arcodrift1_test, na.rm = TRUE)/nrow(arcodrift1_test)*100
# 
# arcodrift1_test_0.001 <- (arcodrift1["p.value",]   < alpha[1])
# arcodrift1_test_0.01 <- (arcodrift1["p.value",]   < alpha[2])
# arcodrift1_test_0.05 <- (arcodrift1["p.value",]   < alpha[3])
# arcodrift1_test_0.10 <- (arcodrift1["p.value",]   < alpha[4])
# arcodrift1_test<- cbind(arcodrift1_test_0.001, arcodrift1_test_0.01, arcodrift1_test_0.05,  arcodrift1_test_0.10)
# test1<- colSums(arcodrift1_test == TRUE)/nrow(arcodrift1_test)*100 #colSums(arcodrift1_test, na.rm = TRUE)/nrow(arcodrift1_test)*100
# 
# arcodrift2_test_0.001 <- (arcodrift2["p.value",]   < alpha[1])
# arcodrift2_test_0.01 <- (arcodrift2["p.value",]   < alpha[2])
# arcodrift2_test_0.05 <- (arcodrift2["p.value",]   < alpha[3])
# arcodrift2_test_0.10 <- (arcodrift2["p.value",]   < alpha[4])
# arcodrift2_test<- cbind(arcodrift2_test_0.001, arcodrift2_test_0.01, arcodrift2_test_0.05,  arcodrift2_test_0.10)
# test2<- colSums(arcodrift2_test == TRUE)/nrow(arcodrift2_test)*100
# 
# TESTS <- rbind(test0,test1,test2)
# colnames(TEST)<- alpha
# rm(test0,test1,test2)
# 
# diff_arco_test_0.001 <- (diff_arco["p.value",]   < alpha[1])
# diff_arco_test_0.01 <- (diff_arco["p.value",]   < alpha[2])
# diff_arco_test_0.05 <- (diff_arco["p.value",]   < alpha[3])
# diff_arco_test_0.10 <- (diff_arco["p.value",]   < alpha[4])
# diff_arco_test<- cbind(diff_arco_test_0.001, diff_arco_test_0.01, diff_arco_test_0.05,  diff_arco_test_0.10)
# difftest0 <- colSums(diff_arco_test == TRUE)/nrow(diff_arco_test)*100
# 
# diff_arcodrift1_test_0.001 <- (diff_arcodrift1["p.value",]   < alpha[1])
# diff_arcodrift1_test_0.01 <- (diff_arcodrift1["p.value",]   < alpha[2])
# diff_arcodrift1_test_0.05 <- (diff_arcodrift1["p.value",]   < alpha[3])
# diff_arcodrift1_test_0.10 <- (diff_arcodrift1["p.value",]   < alpha[4])
# diff_arcodrift1_test<- cbind(diff_arcodrift1_test_0.001, diff_arcodrift1_test_0.01, diff_arcodrift1_test_0.05,  diff_arcodrift1_test_0.10)
# difftest1<- colSums(diff_arcodrift1_test == TRUE)/nrow(diff_arcodrift1_test)*100 
# 
# diff_arcodrift2_test_0.001 <- (diff_arcodrift2["p.value",]   < alpha[1])
# diff_arcodrift2_test_0.01 <- (diff_arcodrift2["p.value",]   < alpha[2])
# diff_arcodrift2_test_0.05 <- (diff_arcodrift2["p.value",]   < alpha[3])
# diff_arcodrift2_test_0.10 <- (diff_arcodrift2["p.value",]   < alpha[4])
# diff_arcodrift2_test<- cbind(diff_arcodrift2_test_0.001, diff_arcodrift2_test_0.01, diff_arcodrift2_test_0.05,  diff_arcodrift2_test_0.10)
# difftest2<- colSums(diff_arcodrift2_test == TRUE)/nrow(diff_arcodrift2_test)*100
# 
# diffTESTS <- rbind(difftest0,difftest1,difftest2)
# colnames(diffTESTS)<- alpha
# rm(difftest0,difftest1,difftest2)
# 
# difftest0 <- percentagesignificant(diff_arco, alpha)
# difftest1 <- percentagesignificant(diff_arcodrift1, alpha)
# difftest2 <- percentagesignificant(diff_arcodrift2, alpha)
# diffTESTS <- rbind(difftest0,difftest1,difftest2)
# colnames(diffTESTS)<- alpha
# rm(test0,test1,test2)