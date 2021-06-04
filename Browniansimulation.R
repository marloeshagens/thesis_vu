### ### ### ### ### SIMULATION EXPERIMENT ### ### ### ### ### 
setwd("/cloud/project")
cat("\014")        
rm(list = ls())

#df <- read_excel("test_production_2.xlsx")
#df<- as.data.frame(df)
#UK <- df[c(35)]
# model <- arima(UK, order = c(1,0,1))
# library(forecast)
#x <- auto.arima(UK)

set.seed(123)
library(ArCo) # contains ArCo package
library(glmnet) #contains LASSO function (Friedman et al., 2010).

#ts.sim <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)
#ts.plot(ts.sim)

######################## BROWNIAN ########################
# https://towardsdatascience.com/monte-carlo-simulation-in-r-with-focus-on-financial-data-ad43e2a4aedf

paths <- 100 # amount of plots/paths we create
T <- 100 # count

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
######################## STANDARD BROWNIAN ########################
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

######################## GEOMETRIC BROWNIAN ########################
GeometricBrownian<-function(paths, count)
{
  #paths<-10
  #count<-5000
  interval<-5/count
  mean<-0.06 #arbitrary; these I can maybe change later, too (and vary the input)
  sigma<-0.3 #arbitrary; these I can maybe change later, too (and vary the input)
  sample<-matrix(0,nrow=(count+1),ncol=paths)
  for(i in 1:paths)
  {
    sample[1,i]<-100
    for(j in 2:(count+1))
    {
      sample[j,i]<-sample[j-1,i]*exp(interval*(mean-((sigma)^2)/2)+((interval)^.5)*rnorm(1)*sigma) #Expression for Geometric Brownian Motion
    }
  }	
  #cat("E[W(2)] = ",mean(sample[2001,]),"\n")
  #cat("E[W(5)] = ",mean(sample[5001,]),"\n")
  matplot(sample,main="Geometric Brownian",xlab="Time",ylab="Path",type="l")
  return(sample)
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
#### data generation ####
x<-Brownian(paths, T)
randomc <- sample(1:paths, 1) #select random path (column) from all created to show
randombrownianx <- x[,randomc]
plot(randombrownianx, type="l")

y<-StandardBrownian(paths, T)
randombrowniany <- y[,randomc] 
plot(randombrowniany, type="l")

z<-GeometricBrownian(paths, T)
randombrownianz <- z[,randomc] 
plot(randombrownianz, type="l")
# randomt0 <- sample(1:24,1) #pick random "treatment" year
randomtreated <- sample(1:paths,1) #pick random "treated" unit

t0=T/2 # I can maybe also use more possibilities?

y<-StandardBrownian(paths, T)
list_y<-list(y) # dataprep arco
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

z<-GeometricBrownian(paths, T)
list_z<-list(z) # dataprep arco
arcoz <- fitArCo(data = list_y,  fn = cv.glmnet, p.fn = predict, 
                 treated.unit = randomtreated, t0 = t0, 
                 boot.cf=TRUE, R=1000 # should I change the bootstrap?
                 #VCOV.type = "nw", prewhitening.kernel = TRUE
)
plot(arcoz, display.fitted=TRUE,confidence.bands = TRUE, alpha = 0.05, main = "Standard Brownian Motion ArCo")
arcoz$p.value  < 0.05
