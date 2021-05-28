# simulation 2 
# https://nl.mathworks.com/help/econ/simulate-trend-stationary-and-difference-stationary-processes.html

t = 1:T
y_stationary <- rnorm(length(t),mean=1,sd=1) # the stationary time series (ts)
y_trend      <- cumsum(rnorm(length(t),mean=1,sd=4))+t/100 # our ts with a trend
# lets normalize each for simplicity
y_stationary<- y_stationary/max(y_stationary) 
y_trend      <- y_trend/max(y_trend) 

adf.test(y_stationary)
kpss.test(y_stationary, null="Trend")
adf.test(y_trend)
kpss.test(y_trend, null="Trend")



#### 
set.seed(123456)
e <- rnorm(500)
rw <- cumsum(e)
trd <- 1:500
rw.wd <- 0.5 * trd + cumsum(e)
dt <- e + 0.5 * trd
ar1 <- arima.sim(model = list(ar = 0.8), n = 500)
plot(dt)