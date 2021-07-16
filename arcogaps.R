# function that extracts data from multiple fitArCo output (listed) and construct gaps 
# note: this function is only useful Y>1 simulations; when Y=1, data is easily extracted

# July 6 2021 
# author: M Hagens

arcogaps <- function(arco_results, diff){
  
  # construct counterfactual
  arco_cf <- rbind(as.data.frame(arco_results["fitted.values",]), as.data.frame(arco_results["cf",]))
  
  # retrieve data and the treated unit (integer) for each simulation
  data <- arco_results["data",]
  treated <- as.data.frame(arco_results["treated.unit",])
  
  # construct Y1 
  if(diff == FALSE){
    Y1 <- matrix(0, T, ncol(arco_results))
    for(i in 1:ncol(arco_results)){
      Y1[,i] <- data[[i]][[1]][,treated[,i]]
    }
  }
  if(diff == TRUE){
    Y1 <- matrix(0, (T-1), ncol(arco_results))
    for(i in 1:ncol(arco_results)){
      Y1[,i] <- data[[i]][[1]][,treated[,i]]
    }
  }
  
  # compute gaps
  arco_gap <- Y1 - arco_cf # each column is one simulation, so nr of columns = Y
}