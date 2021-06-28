# *-----------------------------------------------------------------
# | PROGRAM NAME: ex synthetic controls.R
# | DATE: 3/28/19 
# | CREATED BY: MATT BOGARD 
# | PROJECT FILE: Macintosh HD ▸ ⁨Users⁩ ▸ ⁨amandabogard⁩ ▸ ⁨Google Dive⁩ ▸ ⁨R Training⁩           
# *----------------------------------------------------------------
# | PURPOSE: code based on Synth: An R Package for Synthetic Control Methods in Comparative Case Studies
# | Journal of Statistical Software June 2011, Volume 42, Issue 13
# *----------------------------------------------------------------

# See also: https://cran.r-project.org/web/packages/Synth/Synth.pdf


### Background

# We demonstrate the synthetic control method using data from Abadie and Gardeazabal
# (2003), which studied the economic effects of conflict, using the terrorist conflict in the Basque
# Country as a case study. Abadie and Gardeazabal (2003) used a combination of other Spanish
# regions to construct a synthetic Basque Country resembling many relevant economic characteristics 
# of the Basque Country before the onset of political terrorism in the 1970s.

### Methodology/Estimation

# Synthetic control methods involve the construction of synthetic control units as convex 
# combinations of multiple control units.

#  To construct our synthetic control unit we define a vector of weights W
#  Each W then represents one particular weighted average of control units and therefore one potential synthetic control unit.

# To create the most similar synthetic control unit, the synth() function chooses the vector W∗
# to minimize a distance, kX1 − X0Wk, between X1 and X0W, subject to the weight constraints.

# The V matrix is introduced to allow different weights to the variables in X0 and X1 depending
# on their predictive power on the outcome. An optimal choice of V assigns weights that
# minimize the mean square error of the synthetic control estimator, that is the expectation of
# (Y1 − Y0W∗)'(Y1 − Y0W∗)

# In this procedure a V∗ is chosen among all positive definite and diagonal matrices such that
# the mean squared prediction error (MSPE) of the outcome variable is minimized over some
# set of pre-intervention periods.

### Inference

# Abadie et al. (2010) describe how synthetic control methods facilitate inferential techniques
# akin to permutation tests-so-called placebo studies - iteratively apply the synthetic
# control method by randomly reassigning the intervention in time (i.e., pre-intervention dates)
# or across units (i.e., to control units where the intervention did not occur) to produce a set of
# placebo effects. Subsequently, we can compare the set of placebo effects to the effect that was
# estimated for the time and unit where the intervention actually occurred. This comparison
# is informative about the rarity of the magnitude of the treatment effect that was observed
# for the exposed unit.


library(Synth) # load Synth package

data(basque) # load data

# IMPORTANT NOTES ABOUT THE REQ DATA STRUCTURE

# dataset is organized in standard (long) panel-data format, with variables extending
# across the columns and the rows sorted first by region and then by time-period.10 A name
# (character-string) and number is provided for each region.11 At least one of these two types
# of unit-identifiers is required for Synth

#-------------------------------------
# data prep
#------------------------------------

# The first step is to reorganize the panel dataset into an appropriate format that is suitable
# for the main estimator function synth(). At a minimum, synth() requires as inputs the four
# data matrices X1, X0, Z1, and Z0 that are needed to construct a synthetic control unit. 


dataprep.out <- dataprep(
  foo = basque,
  predictors = c("school.illit", "school.prim", "school.med","school.high", "school.post.high", "invest"),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  special.predictors = list(
    list("gdpcap", 1960:1969 , "mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  controls.identifier = c(2:16, 18),
  time.optimize.ssr = 1960:1969,
  time.plot = 1955:1997)

# dataprep() returns a list object dataprep.out that contains several elements includingX0,X1,Z0,Z1 

print(dataprep.out$X0) # check 
print(dataprep.out$Z1) # Basque GDP per-capita for the pre-intervention period



#------------------------------------------
# data pre-processing/modifications
#-----------------------------------------

# consolidate two fields and calcualte percentage shares
dataprep.out$X1["school.high",] <- dataprep.out$X1["school.high",] +dataprep.out$X1["school.post.high",]
dataprep.out$X1 <- as.matrix(dataprep.out$X1[-which(rownames(dataprep.out$X1) == "school.post.high"),])
dataprep.out$X0["school.high",] <- dataprep.out$X0["school.high",] +dataprep.out$X0["school.post.high",]
dataprep.out$X0 <- dataprep.out$X0[-which(rownames(dataprep.out$X0) == "school.post.high"),]
lowest <- which(rownames(dataprep.out$X0) == "school.illit")
highest <- which(rownames(dataprep.out$X0) == "school.high")
dataprep.out$X1[lowest:highest,] <- (100 * dataprep.out$X1[lowest:highest,]) /sum(dataprep.out$X1[lowest:highest,])


#-------------------------------------
# analysis using synth
#-------------------------------------

# The synth() command searches for the W∗ vector of weights that identifies the synthetic control 
# for the Basque region by solving the nested optimization problem 

synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")

# calculate GPD difference in trend between the Basque region and its synthetic control
gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth.out$solution.w)
print(gaps[1:3, 1]) # check

# pre built tables from synth objects
synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
names(synth.tables) # check

# comparing pre-treatment predictor values for the treated unit, the synthetic control unit, and all the units in the sample

synth.tables$tab.pred[1:5, ] # check balance across treated and control for pre-period predictors

# view control unit weights
synth.tables$tab.w[8:14, ]

# plot treatment vs control outcomes for pre and post periods
path.plot(synth.res = synth.out, dataprep.res = dataprep.out,
          Ylab = "real per-capita GDP (1986 USD, thousand)", Xlab = "year",
          Ylim = c(0, 12), Legend = c("Basque country","synthetic Basque country"), 
          Legend.position = "bottomright")

# gap plot
gaps.plot(synth.res = synth.out, dataprep.res = dataprep.out,
          Ylab = "gap in real per-capita GDP (1986 USD, thousand)", Xlab = "year",
          Ylim = c(-1.5, 1.5), Main = NA)

#--------------------------------------
# inference via placebo tests
#-------------------------------------

# These tests involve applying the synthetic control method after reassigning the intervention in the
# data to units and periods where the intervention did not occur

store <- matrix(NA,length(1955:1997),17)
colnames(store) <- unique(basque$regionname)[-1]

# run placebo test
for(iter in 2:18)
{
  dataprep.out <-
    dataprep(foo = basque,
             predictors = c("school.illit" , "school.prim" , "school.med" ,
                            "school.high" , "school.post.high" , "invest") ,
             predictors.op = "mean" ,
             time.predictors.prior = 1964:1969 ,
             special.predictors = list(
               list("gdpcap" , 1960:1969 , "mean"),
               list("sec.agriculture" ,      seq(1961,1969,2), "mean"),
               list("sec.energy" ,           seq(1961,1969,2), "mean"),
               list("sec.industry" ,         seq(1961,1969,2), "mean"),
               list("sec.construction" ,     seq(1961,1969,2), "mean"),
               list("sec.services.venta" ,   seq(1961,1969,2), "mean"),
               list("sec.services.nonventa" ,seq(1961,1969,2), "mean"),
               list("popdens", 1969, "mean")
             ),
             dependent = "gdpcap",
             unit.variable = "regionno",
             unit.names.variable = "regionname",
             time.variable = "year",
             treatment.identifier = iter,
             controls.identifier = c(2:18)[-iter+1],
             time.optimize.ssr = 1960:1969,
             time.plot = 1955:1997
    )
  
  
  dataprep.out$X1["school.high",] <-dataprep.out$X1["school.high",] + dataprep.out$X1["school.post.high",]
  dataprep.out$X1 <- as.matrix(dataprep.out$X1[-which(rownames(dataprep.out$X1)=="school.post.high"),])
  dataprep.out$X0["school.high",] <-dataprep.out$X0["school.high",] + dataprep.out$X0["school.post.high",]
  dataprep.out$X0 <-dataprep.out$X0[-which(rownames(dataprep.out$X0)=="school.post.high"),]
  
  lowest  <- which(rownames(dataprep.out$X0)=="school.illit")
  highest <- which(rownames(dataprep.out$X0)=="school.high")
  
  dataprep.out$X1[lowest:highest,] <-(100*dataprep.out$X1[lowest:highest,]) /sum(dataprep.out$X1[lowest:highest,])
  dataprep.out$X0[lowest:highest,] <-100*scale(dataprep.out$X0[lowest:highest,],center=FALSE,scale=colSums(dataprep.out$X0[lowest:highest,])
  )
  
  # run synth
  synth.out <- synth(
    data.prep.obj = dataprep.out,
    method = "BFGS"
  )
  
  # store gaps
  store[,iter-1] <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth.out$solution.w)
}

# now do figure
data <- store
rownames(data) <- 1955:1997

# Set bounds in gaps data
gap.start     <- 1
gap.end       <- nrow(data)
years         <- 1955:1997
gap.end.pre  <- which(rownames(data)=="1969")

#  MSPE Pre-Treatment
mse        <- apply(data[ gap.start:gap.end.pre,]^2,2,mean)
basque.mse <- as.numeric(mse[16])


# Exclude states with 5 times higher MSPE than basque
data <- data[,mse<5*basque.mse]
Cex.set <- .75

# Plot
plot(years,data[gap.start:gap.end,which(colnames(data)=="Basque Country (Pais Vasco)")],
     ylim=c(-2,2),xlab="year",
     xlim=c(1955,1997),ylab="gap in real per-capita GDP (1986 USD, thousand)",
     type="l",lwd=2,col="black",
     xaxs="i",yaxs="i")

# Add lines for control states
for (i in 1:ncol(data)) { lines(years,data[gap.start:gap.end,i],col="gray") }

## Add Basque Line
lines(years,data[gap.start:gap.end,which(colnames(data)=="Basque Country (Pais Vasco)")],lwd=2,col="black")

# Add grid
abline(v=1970,lty="dotted",lwd=2)
abline(h=0,lty="dashed",lwd=2)
legend("bottomright",legend=c("Basque country","control regions"),
       lty=c(1,1),col=c("black","gray"),lwd=c(2,1),cex=.8)
arrows(1967,-1.5,1968.5,-1.5,col="black",length=.1)
text(1961.5,-1.5,"Terrorism Onset",cex=Cex.set)
abline(v=1955)
abline(v=1997)
abline(h=-2)
abline(h=2)