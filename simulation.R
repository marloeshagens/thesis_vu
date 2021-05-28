### ### ### ### ### SIMULATION EXPERIMENT ### ### ### ### ### 
df <- read_excel("test_production_2.xlsx")
df<- as.data.frame(df)
UK <- df[c(35)]

# model <- arima(UK, order = c(1,0,1))
library(forecast)
x <- auto.arima(UK)

set.seed(123)
x.sim <-  arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)

ts.sim <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)
ts.plot(ts.sim)


######################## BROWNIAN ########################
# https://towardsdatascience.com/monte-carlo-simulation-in-r-with-focus-on-financial-data-ad43e2a4aedf

paths <- 100 # amount of plots/paths we create
T <- 29 # count

Brownian<-function(paths, count) # This is a function to generate Brownian with drift 0.06 and diffusion 0.3
{
  #paths<-100
  #count<-5000
  interval<-5/count
  sample<-matrix(0,nrow=(count+1),ncol=paths)
  for(i in 1:paths)
  {
    sample[1,i]<-5
    for(j in 2:(count+1))
    {
      sample[j,i]<-sample[j-1,i]+interval*0.06+((interval)^.5)*rnorm(1)*0.3
    }
  }	
  #cat("E[W(2)] = ",mean(sample[2001,]),"\n")
  #cat("E[W(5)] = ",mean(sample[5001,]),"\n")
  matplot(sample,main="Brownian",xlab="Time",ylab="Path",type="l")
  return(sample)
}
StandardBrownian<-function(paths, count) # This is a function to generate Standard Browninan with drift 0 and diffusion 1
{
  #paths<-10
  #count<-5000
  interval<-5/count
  sample<-matrix(0,nrow=(count+1),ncol=paths)
  for(i in 1:paths)
  {
    sample[1,i]<-0
    for(j in 2:(count+1))
    {
      sample[j,i]<-sample[j-1,i]+((interval)^.5)*rnorm(1)
    }
  }	
  #cat("E[W(2)] = ",mean(sample[2001,]),"\n")
  #cat("E[W(5)] = ",mean(sample[5001,]),"\n")
  matplot(sample,main="Standard Brownian",xlab="Time",ylab="Path",type="l")
  return(sample)
}

x<-Brownian(paths, T)
randomc <- sample(1:paths, 1) #select random path from all created
randombrownianx <- x[,randomc]
plot(randombrownianx, type="l")

y<-StandardBrownian(paths, T)
randombrowniany <- y[,randomc] #randomc=15
plot(randombrowniany, type="l")

# randomt0 <- sample(1:24,1) #pick random "treatment" year
randomtreated <- sample(1:paths,1) #pick random "treated" unit

list_y<-list(y)
arcoy <- fitArCo(data = list_y,  fn = cv.glmnet, p.fn = predict, 
                treated.unit = randomtreated, t0 = t0, 
                boot.cf=TRUE,R=200,l=3 
                #VCOV.type = "nw", prewhitening.kernel = TRUE
)
plot(arcoy, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Standard Brownian Motion ArCo")
arcoy$p.value  < 0.05

adf.test(y[,randomtreated])
p.valuekpss <- kpss.test(y[,randomtreated], null="Trend")$p.value
p.valuekpss <0.05


GeometricBrownian<-function(paths, count)
{
  #paths<-10
  #count<-5000
  interval<-5/count
  mean<-0.06
  sigma<-0.3
  sample<-matrix(0,nrow=(count+1),ncol=paths)
  for(i in 1:paths)
  {
    sample[1,i]<-100
    for(j in 2:(count+1))
    {
      sample[j,i]<-sample[j-1,i]*exp(interval*(mean-((sigma)^2)/2)+((interval)^.5)*rnorm(1)*sigma) #Expression for Geometric Brownian Motion
    }
  }	
  cat("E[W(2)] = ",mean(sample[2001,]),"\n")
  cat("E[W(5)] = ",mean(sample[5001,]),"\n")
  matplot(sample,main="Geometric Brownian",xlab="Time",ylab="Path",type="l")
}
GeometricBrownian()