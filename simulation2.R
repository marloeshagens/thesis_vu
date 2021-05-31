# simulation 2 
# https://nl.mathworks.com/help/econ/simulate-trend-stationary-and-difference-stationary-processes.html

set.seed(12345)
T=29
t = 1:T
y_stationary <- rnorm(length(t),mean=1,sd=1) # the stationary time series (ts)
y_trend      <- cumsum(rnorm(length(t),mean=-0.5,sd=2))+0.5*t/100 # our ts with a trend
#repeat above equation multiple times
y<- ts(replicate(n = 5, cumsum(rnorm(length(t),mean=-0.5,sd=2))+0.5*t/100))

# lets normalize each for simplicity
y_stationary<- y_stationary/max(y_stationary) 
y_trend      <- y_trend/max(y_trend) 

adf.test(y_stationary)
kpss.test(y_stationary, null="Trend")
adf.test(y_trend)
kpss.test(y_trend, null="Trend")

plot(y_trend)
plot(y_stationary)

par(mfcol=c(1,1))
plot(y)
matplot(y, 
        type = "l", 
        col = c("steelblue", "darkgreen", "darkred", "orange","yellow"), 
        lty = 1, 
        lwd = 2,
        main = "ARIMA (1,1,0), ar=0.7, mean =-0.5",
        xlab = "Time",
        ylab = "Value")

plot(t,y_stationary,
     type='l',col='red',
     xlab = "time (t)",
     ylab = "Y(t)",
     main = "Stationary signal")
acf(y_stationary,lag.max = length(y_stationary),
    xlab = "lag #", ylab = 'ACF',main=' ')
# the trend signal and ACF
plot(t,y_trend,
     type='l',col='red',
     xlab = "time (t)",
     ylab = "Y(t)",
     main = "Trend signal")
acf(y_trend,lag.max = length(y_trend),
    xlab = "lag #", ylab = 'ACF', main=' ')

#### 
# set.seed(123456)
# e <- rnorm(500)
# rw <- cumsum(e)
# trd <- 1:500
# rw.wd <- 0.5 * trd + cumsum(e)
# dt <- e + 0.5 * trd
# ar1 <- arima.sim(model = list(ar = 0.8), n = 500)
# plot(dt)

###
t = 0:300
y_stationary <- rnorm(length(t),mean=1,sd=1) # the stationary time series (ts)
y_trend      <- cumsum(-rnorm(length(t),mean=-1,sd=1))+0.5*t/100 # our ts with a trend
# lets normalize each for simplicity
y_stationary<- y_stationary/max(y_stationary) 
y_trend      <- y_trend/max(y_trend) 
plot.new()
frame()
par(mfcol=c(2,2))
# the stationary signal and ACF
plot(t,y_stationary,
     type='l',col='red',
     xlab = "time (t)",
     ylab = "Y(t)",
     main = "Stationary signal")
acf(y_stationary,lag.max = length(y_stationary),
    xlab = "lag #", ylab = 'ACF',main=' ')
# the trend signal and ACF
plot(t,y_trend,
     type='l',col='red',
     xlab = "time (t)",
     ylab = "Y(t)",
     main = "Trend signal")
acf(y_trend,lag.max = length(y_trend),
    xlab = "lag #", ylab = 'ACF', main=' ')

####
# simulate and plot ARIMA (1,1,0)
set.seed(12345)
trend_stat_ts <- ts(replicate(n = 5, 
                     arima.sim(model = list(order = c(1, 1, 0), ar = 0.7),
                               n = 29,
                               mean = -0.5)))
par(mfcol=c(1,1))
matplot(trendstatts, 
        type = "l", 
        col = c("steelblue", "darkgreen", "darkred", "orange"), 
        lty = 1, 
        lwd = 2,
        main = "ARIMA (1,1,0), ar=0.7, mean =-0.5",
        xlab = "Time",
        ylab = "Value")
