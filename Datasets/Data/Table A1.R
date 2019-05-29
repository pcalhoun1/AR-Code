#rm(list=ls(all=TRUE))

library(dplyr)

load(file="binData.RData")
load(file="contData.RData")


### Binary Outcomes ###
binTable <- data.frame(dataset = names(binData), nObs = unlist(lapply(binData, nrow)), nTotalVars = unlist(lapply(binData, ncol))-1, stringsAsFactors = FALSE)

# Not all variables are used.  Must go through each dataset and change nTotalVars based on the formula used to grow the RF
# For simplicity, only show the changes, but checked all of them
str(binData$breastCancer)
# Note: ID not used
binTable[binTable$dataset == "breastCancer", 3] <- binTable[binTable$dataset == "breastCancer", 3] - 1

# Note: V2 only has one level, not used in RF 
str(binData$ionosphere)
binTable[binTable$dataset == "ionosphere", 3] <- binTable[binTable$dataset == "ionosphere", 3] - 1


### Continuous Outcome ###
contTable <- data.frame(dataset = names(contData), nObs = unlist(lapply(contData, nrow)), nTotalVars = unlist(lapply(contData, ncol))-1, stringsAsFactors = FALSE)
contTable$multiplier <- "1"
contTable$multiplier[contTable$dataset == "servo"] <- "10"
contTable$multiplier[contTable$dataset == "friedman2"] <- "10^4"
contTable$multiplier[contTable$dataset == "friedman3"] <- "10^(-2)"
contTable$multiplier[contTable$dataset == "ailerons"] <- "10^(-8)"
contTable$multiplier[contTable$dataset == "elevators"] <- "10^(-6)"
contTable$multiplier[contTable$dataset == "imports85"] <- "10^6"
contTable$multiplier[contTable$dataset == "airquality"] <- "10^2"

str(contData$housing)

# Table A1
bind_rows(binTable, contTable)

