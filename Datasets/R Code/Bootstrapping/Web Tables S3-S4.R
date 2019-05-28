#rm(list=ls(all=TRUE))

library(dplyr)
library(reshape2)

##### BOOTSTRAP #####

fileList <- list.files("minpvalue/Results")

binResults <- NULL
contResults <- NULL
for (file in fileList){
  tempDat <-  read.csv(paste0("minpvalue/Results/", file))
  if (names(tempDat)[5] == "misclass"){
    binResults <- rbind(binResults, tempDat)
  } else if (names(tempDat)[5] == "mse") {
    contResults <- rbind(contResults, tempDat)
  } else { stop("Unexpected accuracy metric") }
}
dim(binResults)
dim(contResults)

# Add in minpvalue = 0.05 #
fileListReg <- list.files("Individual Dataset/Results")
for (file in fileListReg){
  tempDat <-  read.csv(paste0("Individual Dataset/Results/", file))
  tempDat$minpvalue <- 0.05
  if (names(tempDat)[4] == "misclass"){
    binResults <- rbind(binResults, tempDat[tempDat$method=="AR", ])
  } else if (names(tempDat)[4] == "mse") {
    contResults <- rbind(contResults, tempDat[tempDat$method=="AR", ])
  } else { stop("Unexpected accuracy metric") }
}
dim(binResults)
dim(contResults)



### Binary Datasets ###
head(binResults)

### AUC ###
binAUC <- aggregate(list(binResults[["auc"]]), list(binResults[["dataset"]], binResults[["method"]], binResults[["minpvalue"]]), FUN=function(x){c(mean=round(mean(x),4), sd=round(sd(x),5))})
binAUC <- data.frame(dataset=binAUC[,1], method=binAUC[,2], minpvalue=binAUC[,3], aucMean=binAUC[,4][,1], aucSD = binAUC[,4][,2])

### Table S3 ###
(aucMean <- dcast(binAUC, dataset ~ minpvalue + method, value.var="aucMean")[c(1,7,4,3,5,10,9,6,8,2), ])  # Mean
#To compare AUC, calculate 1-AUC.  The difference between 0.98 and 0.99 is further than 0.01 and 0.02
table(factor(unlist(apply(aucMean[2:ncol(aucMean)], 1, function(x) { names(x)[x == max(x)] })), levels=c("0.01_AR","0.05_AR","0.1_AR","0.25_AR","0.5_AR","0.9_AR")))  #Number of times method was most accurate
bestResult <- apply(1 - aucMean[2:ncol(aucMean)], 1, min)
round(colMeans(((1 - aucMean[2:ncol(aucMean)]) - bestResult) / bestResult)*100, 1)

### SD of AUC (data not shown) ###
dcast(binAUC, dataset ~ minpvalue + method, value.var="aucSD")[c(1,7,4,3,5,10,9,6,8,2), ]  # SD



### Continuous Datasets ###

# Create multiplier
contResults$multiplier <- 1
contResults$multiplier[contResults$dataset == "Servo"] <- 10
contResults$multiplier[contResults$dataset == "Friedman2"] <- 10^4
contResults$multiplier[contResults$dataset == "Friedman3"] <- 10^(-2)
contResults$multiplier[contResults$dataset == "Ailerons"] <- 10^(-8)
contResults$multiplier[contResults$dataset == "Elevators"] <- 10^(-6)
contResults$multiplier[contResults$dataset == "Imports85"] <- 10^6
contResults$multiplier[contResults$dataset == "Airquality"] <- 10^2

contResults$mseAdj <- contResults$mse * (contResults$multiplier)^-1


### MSE ###
contMSE <- aggregate(list(contResults[["mseAdj"]]), list(contResults[["dataset"]], contResults[["method"]], contResults[["minpvalue"]]), FUN=function(x){c(mean=round(mean(x),3), sd=round(sd(x),3))})
contMSE <- data.frame(dataset=contMSE[,1], method=contMSE[,2], minpvalue=contMSE[,3], mseMean=contMSE[,4][,1], mseSD = contMSE[,4][,2])

### Table S4 ###
(mseMean <- dcast(contMSE, dataset ~ minpvalue + method, value.var="mseMean")[c(8,10,1,5,6,7,2,4,9,3), ])  # Mean
table(factor(unlist(apply(mseMean[2:ncol(mseMean)], 1, function(x) { names(x)[x == min(x)] })), levels=c("0.01_AR","0.05_AR","0.1_AR","0.25_AR","0.5_AR","0.9_AR")))  #Number of times method was most accurate
bestResult <- apply(mseMean[2:ncol(mseMean)], 1, min)
round(colMeans(((mseMean[2:ncol(mseMean)]) - bestResult) / bestResult)*100, 1)

### SD of MSE (data not shown) ###
dcast(contMSE, dataset ~ minpvalue + method, value.var="mseSD")[c(8,10,1,5,6,7,2,4,9,3), ]  # SD


