##############################################################################
#### function for the nonstationary time series counterfactual experiment ####
##############################################################################
# July 6 2021 
# author: M Hagens

### summary ###
# this function creates a N*T matrix of random walk data,
# puts in into a counterfactual computing model
# and gives counterfactual results 

### input options ### 
# model: either "ArCo" or "SC" #
# T:  the length of each time series #
# N: the number of units #
# delta: the (possible) drift in the random walk # 
# variance: specifies the random walk's error term distribution #
# R: if model = "ArCo", R specifies the number of bootstrap replications in the ArCo #

### the output of this function ### 
# for model="ArCo": the complete fitArCo results in a list, incl data, counterfactual, delta's 
# for model="SC": the complete data on which the counterfactual is build, incl the treated unit, the counterfactual, and gaps between counterfactual and data

nonstationary_ts_counterfactual <- function(model, T, N, y0, delta, variance, R, diff){
  
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
  
  # model selection # 
  if(model == "ArCo"){ # put dataframe in ArCo format, run and return ArCo
    
    list<-list(RW)
    
    arco <- fitArCo(
      data = list,  
      fn = cv.glmnet, 
      p.fn = predict, 
      treated.unit = randomtreated, 
      t0 = t0, 
      boot.cf = TRUE, 
      R = R,
      VCOV.type = "iid"
    )
    
    plot(arco, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Random Walk ArCo")
    
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
    dataprep.out<- dataprep(
      foo = dfLong,
      predictors = "X", # in this case, we predict the dependent variable on its own values.
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
    
    # run synth # 
    synth.out <- synth(dataprep.out, method = "BFGS")
    
    # plot treatment vs control outcomes for pre and post periods #
    path.plot(synth.res = synth.out, dataprep.res = dataprep.out,
              Ylab = "Value", Xlab = "Time",
              Legend = c("Random Walk","Synthetic Control"), 
              Legend.position = "bottomright",
              Main = "Simulated RW data and its SC")
    
    # gap plot #
    gaps.plot(synth.res = synth.out, dataprep.res = dataprep.out,
              Ylab = "Value", Xlab = "Time",
              Main = "Gap between RW and SC")
    
    # make output #
    Y0plot<- dataprep.out$Y0plot
    Y1 <- dataprep.out$Y1plot
    counterfactual <- dataprep.out$Y0plot%*%synth.out$solution.w
    gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot%*%synth.out$solution.w) #difference in trend between simulated data and its synthetic control
    
    SC_output<- (cbind.data.frame(Y0plot, Y1, counterfactual, gaps))
    colnames(SC_output)[20:22]<- c("Y1", "counterfactual","gaps")
    
    return(SC_output)
    
    # you can change the return if you are interested in more general results of the SC,
    # namely the solutions for predictor and control weights and loss from optimization (MSPE) and the 
    # results from the optimx() minimization, which could be used for diagnostics. 
    # return(synth.out) 
    
    ls()
  }
  ###########
  ### END ### 
  ###########
}