### SYNTH CONTROL on CO2 EMISSION DATA 

#load data
df <- read.csv("historical_emissions.csv", header = TRUE)
df <- df[-c(2:5)]
for (col in 2:ncol(df)){
  colnames(df)[col] <-  sub("X","", colnames(df)[col])
}

# put DF in SC format

# the melt() function can reshape a long data frame to wide data frame, 
# which is needed as input for SC
dfLong <- melt(df, id.vars="Country", 
               variable.name = "Year", value.name="MtCo2e")
# head(dfLong)

# create a numeric unit variable to indicate which country is which "unit.num"
N <- nrow(df) # N = nr of countries, incl treated unit.
T <- nrow(tdf) # T = nr of Years
dfLong$unit.num <- rep(c(1:N),times=T)

# create matrices from panel data that provide inputs for synth()
dataprep.out<-
  dataprep(
    foo = dfLong,
    predictors = "MtCo2e",
    predictors.op = "mean", #mean is the default
    dependent = "MtCo2e",
    unit.variable = "unit.num",  # this has to be a numeric variable
    time.variable = "Year",
          # not sure what the below command is for.
          # special.predictors = list( 
           #  list("Y", 1991, "mean"),
            # list("Y", 1985, "mean"),
            # list("Y", 1980, "mean")
          # ),
    treatment.identifier = 29, #united kingdom
    controls.identifier = c(1:28), # you have to specify at least two control units
    time.predictors.prior = c(1990:2012),
    time.optimize.ssr = c(1990:2013), 
    unit.names.variable = "Country", 
    time.plot = 1990:2018
  )

# run synth
# identifies weights that create the best possible synthetic control unit for the treated.
synth.out <- synth(dataprep.out)

## summarizing the results by accessing the output from synth.out directly
round(synth.out$solution.w,3) # contains the unit weights
synth.out$solution.v  ## contains the predictor weights. 
# MH: here I only have 1 predictor so weight = 1

## the output from synth opt 
## can be flexibly combined with 
## the output from dataprep to 
## compute other quantities of interest
## for example, the period by period 
## discrepancies between the 
## treated unit and its synthetic control unit
## can be computed by typing
gaps<- dataprep.out$Y1plot-(
  dataprep.out$Y0plot%*%synth.out$solution.w
) ; gaps

# Get result tables - V and W weights plus balance btw. treated and SC
synth.tables <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res = synth.out)
print(synth.tables)

## summary plots for outcome trajectories of the treated and the SC unit

## plot in levels (path; treated and synthetic)
path.plot(dataprep.res = dataprep.out, 
          synth.res = synth.out,
          Ylab = c("CO2 emissions (MtCO2e)"),
          Xlab = c("Year"), 
          #Ylim = c(0,13), 
          Legend = c("UK","Synthetic UK"),
) 
## plot the gaps (treated - synthetic)
gaps.plot(dataprep.res = dataprep.out,
          synth.res = synth.out,
          Ylab = c("Gap in CO2 emissions (MtCO2e)"),
          Xlab = c("Year"))