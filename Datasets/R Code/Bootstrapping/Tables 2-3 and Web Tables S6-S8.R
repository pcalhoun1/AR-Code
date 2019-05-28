#rm(list=ls(all=TRUE))

library(dplyr)
library(reshape2)

##### BOOTSTRAP #####

fileList <- list.files("Individual Dataset/Results")

binResults <- NULL
contResults <- NULL
for (file in fileList){
  tempDat <-  read.csv(paste0("Individual Dataset/Results/", file))
  if (names(tempDat)[4] == "misclass"){
    binResults <- rbind(binResults, tempDat)
  } else if (names(tempDat)[4] == "mse") {
    contResults <- rbind(contResults, tempDat)
  } else { stop("Unexpected accuracy metric") }
}
dim(binResults)
dim(contResults)

### Binary Datasets ###
head(binResults)


### AUC ###
binAUC <- aggregate(list(binResults[["auc"]]), list(binResults[["dataset"]], binResults[["method"]]), FUN=function(x){c(mean=round(mean(x),4), sd=round(sd(x),5))})
binAUC <- data.frame(dataset=binAUC[,1], method=binAUC[,2], aucMean=binAUC[,3][,1], aucSD = binAUC[,3][,2])

### Table 2 ###
(aucMean <- dcast(binAUC, dataset ~ method, value.var="aucMean")[c(1,7,4,3,5,10,9,6,8,2), c(1,4,5,3,2)])  # Mean
#To compare AUC, calculate 1-AUC.  The difference between 0.98 and 0.99 is further than 0.01 and 0.02
table(factor(unlist(apply(aucMean[c("RF","SSS","ER","AR")], 1, function(x) { names(x)[x == max(x)] })), levels=c("RF","SSS","ER","AR")))  #Number of times method was most accurate
bestResult <- apply(1 - aucMean[c("RF","SSS","ER","AR")], 1, min)
round(colMeans(((1 - aucMean[c("RF","SSS","ER","AR")]) - bestResult) / bestResult)*100, 1)

### Table S6 ###
dcast(binAUC, dataset ~ method, value.var="aucSD")[c(1,7,4,3,5,10,9,6,8,2) , c(1,4,5,3,2)]  # SD



### Misclassification ###
binMisclass <- aggregate(list(binResults[["misclass"]]), list(binResults[["dataset"]], binResults[["method"]]), FUN=function(x){c(mean=round(mean(x)*100,2), sd=round(sd(x),6))})
binMisclass <- data.frame(dataset=binMisclass[,1], method=binMisclass[,2], misclassMean=binMisclass[,3][,1], misclassSD = binMisclass[,3][,2])

### Table S8 ###
(misclassMean <- dcast(binMisclass, dataset ~ method, value.var="misclassMean")[c(1,7,4,3,5,10,9,6,8,2) , c(1,4,5,3,2)])  # Mean
table(factor(unlist(apply(misclassMean[c("RF","SSS","ER","AR")], 1, function(x) { names(x)[x == min(x)] })), levels=c("RF","SSS","ER","AR")))  #Number of times method was most accurate
bestResult <- apply(misclassMean[c("RF","SSS","ER","AR")], 1, min)
round(colMeans((misclassMean[c("RF","SSS","ER","AR")]-bestResult)/bestResult)*100, 1)

### SD of Misclassification (data not shown) ###
dcast(binMisclass, dataset ~ method, value.var="misclassSD")[c(1,7,4,3,5,10,9,6,8,2) , c(1,4,5,3,2)]  # SD



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
contMSE <- aggregate(list(contResults[["mseAdj"]]), list(contResults[["dataset"]], contResults[["method"]]), FUN=function(x){c(mean=round(mean(x),3), sd=round(sd(x),3))})
contMSE <- data.frame(dataset=contMSE[,1], method=contMSE[,2], mseMean=contMSE[,3][,1], mseSD = contMSE[,3][,2])

### Table 3 ###
(mseMean <- dcast(contMSE, dataset ~ method, value.var="mseMean")[c(8,10,1,5,6,7,2,4,9,3), c(1,4,5,3,2)])  # Mean
table(factor(unlist(apply(mseMean[c("RF","SSS","ER","AR")], 1, function(x) { names(x)[x == min(x)] })), levels=c("RF","SSS","ER","AR")))  #Number of times method was most accurate
bestResult <- apply(mseMean[c("RF","SSS","ER","AR")], 1, min)
round(colMeans((mseMean[c("RF","SSS","ER","AR")]-bestResult)/bestResult)*100, 1)

### Table S7 ###
dcast(contMSE, dataset ~ method, value.var="mseSD")[c(8,10,1,5,6,7,2,4,9,3) , c(1,4,5,3,2)]  # SD


