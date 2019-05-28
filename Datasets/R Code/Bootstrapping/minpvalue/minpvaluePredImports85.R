#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")


### Imports85 ###

nsim=50

ntrees=100

# imports85 #
imports85<-contData$imports85
form<-as.formula("Response ~ symboling + make + fuelType + aspiration + numOfDoors + bodyStyle +
                 driveWheels + engineLocation + wheelBase + length + width + height + curbWeight +
                 engineType + numOfCylinders + engineSize + fuelSystem + bore + stroke +
                 compressionRatio + horsepower + peakRpm + cityMpg + highwayMpg")

mseAR_imports85 <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_imports85<-growRF_Parallel(ntrees=ntrees, formula=form, data=imports85, search="ar", method="anova", split="MSE",
                                  mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_imports85[index] <- mean((predictRF(ar_imports85,imports85,checkCases=TRUE)-imports85$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Imports85", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_imports85)

#write.table(outData, file = "Results/minpvaluePredImports85.csv", sep = ",", row.names=FALSE)





