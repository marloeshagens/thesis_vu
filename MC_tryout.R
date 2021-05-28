library(MonteCarlo)
library(forecast)
library(Metrics)

# parameter grids
n <- 10 # length of time series
lb <- seq(n-2) + 1 # vector of block sizes
phi <- 0.6 # autoregressive parameter
reps <- 3 # monte carlo replications

# simulation function  
bootstrap1 <- function(n, lb, phi) {
  
  #### simulate ####
  ts <- arima.sim(n, model = list(ar = phi, order = c(1, 0, 0)), sd = 1)
  
  #### divide ####
  m <- ceiling(n / lb) # number of blocks
  blk <- split(ts, rep(1:m, each = lb, length.out = n)) # divide into blocks
  #### resample ####
  res <- sample(blk, replace = TRUE, 10)        # resamples the blocks
  res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series
  #### train, forecast ####
  train <- head(res.unlist, round(length(res.unlist) - 10)) # train set
  test <- tail(res.unlist, length(res.unlist) - length(train)) # test set
  nfuture <- forecast(train, # forecast
                      model = auto.arima(train), 
                      lambda = 0, biasadj = TRUE, h = length(test))$mean    
  ### metric ####
  RMSE <- rmse(test, nfuture) # return RMSE
  return(
    list("RMSE" = RMSE)
  )
}

param_list = list("n" = n, "lb" = lb, "phi" = phi)


set.seed(123)
MC_result <- MonteCarlo(func = bootstrap1, 
                        nrep = reps,
                        ncpus = parallel::detectCores() - 1,
                        param_list = param_list,
                        export_also = list(
                          "packages" = c("forecast", "Metrics")
                        ),
                        raw = T)
