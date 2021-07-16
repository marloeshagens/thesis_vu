# function that computes the percentages of arco outputs that are significant
# July 6 2021 
# author: M Hagens

percentagesignificant <- function(arcooutput, alpha){
  test <- matrix(NA, length(arcooutput["p.value",]), length(alpha))
  for (i in 1:length(alpha)) {
    test[,i] <- (arcooutput["p.value",] < alpha[i])
  }
  percentage_true <- colSums(test == TRUE)/nrow(test)
  return(percentage_true)
}