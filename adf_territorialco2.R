# Territorial emissions in MtCO₂
#Territorial	"Cite as: Friedlingstein et al. (2020).
#Friedlingstein et al., 2020 : The Global Carbon Budget 2020, Earth System Science Data. Available at: <a href=""https://doi.org/10.5194/essd-12-3269-2020"" target=""_blank"">https://doi.org/10.5194/essd-12-3269-2020</a>
#The use of data is conditional on citing the original data sources. 

#Data for years 2018 and 2019 are preliminary."

setwd("/cloud/project")
cat("\014")        
rm(list = ls())
################################################################################  
############################### 0. LOAD PACKAGES ############################### 
################################################################################ 

library(xlsx)
library(ggplot2)
library(moments)
library(dplyr)

# time series
library(tseries)

# Synthetic Control method
install.packages("Synth")
library(Synth)

# Artificial Counterfactual method
install.packages("ArCo")
library(ArCo)

################################################################################  

emissions_data <- read.xlsx("carbonatlas_ts_gb_eu.xlsx", 1)
emissions_data <- subset(emissions_data, select=-c(Year))
emissions_data <- list(emissions_data)
arco <- fitArCo(data = emissions_data, fn = fn, p.fn = p.fn, treated.unit = 3 , t0 = 33)
plot(arco)

#EU27
df_EU27 <- emissions_data %>% select(1,3)

plot(ts(df_EU27$EU27, start = 1960 , end = 2019),
     ylab = "Territorial emissions (MtCO₂)", xlab = "Year",
     type = "l", 
     xlim=c(1960,2019))

adf.test(df_EU27$EU27) # the data set has a unit root.

slimmed_df <- df_EU27[-c(1:30), ]
plot(ts(df_EU27$EU27, start = 1990 , end = 2019),
     ylab = "Territorial emissions (MtCO₂)", xlab = "Year",
     type = "l", 
     xlim=c(1990,2019))
adf.test(slimmed_df$EU27) # the data set has a unit root.

#GB
df_GB <- emissions_data %>% select(1,4)

plot(ts(df_GB$United.Kingdom, start = 1960 , end = 2019),
     ylab = "Territorial emissions (MtCO₂)", xlab = "Year",
     type = "l", 
     xlim=c(1960,2019))

adf.test(df_GB$United.Kingdom) # the data set has a unit root.

#EU excl GB
emissions_data$EU_exclGB <- (emissions_data$EU27 - emissions_data$United.Kingdom)
df_EU_exclGB <- emissions_data %>% select(1,5)

plot(ts(emissions_data$EU_exclGB, start = 1960 , end = 2019),
     ylab = "Territorial emissions (MtCO₂)", xlab = "Year",
     type = "l", 
     xlim=c(1960,2019))

adf.test(emissions_data$EU_exclGB) # the data set has a unit root.
