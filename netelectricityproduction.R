y <- read.csv("netelectricityproduction.tsv", sep = "\t")

df<- as.data.frame(y)
library(dplyr)
dfchar %>% separate(df[,1], c("siec","unit","geo"), sep = "([,])")
#>      A    B
#> 1 <NA> <NA>
#> 2    x    y
#> 3    x    z
#> 4    y    z


data.frame(do.call("rbind", strsplit(as.character(df[,1]), ",", fixed = TRUE)))
