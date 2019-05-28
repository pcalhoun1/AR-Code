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

mseRF_imports85 <- rep(NA, nsim); mseSSS_imports85 <- rep(NA, nsim); mseER_imports85 <- rep(NA, nsim); mseAR_imports85 <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_imports85<-growRF_Parallel(ntrees=ntrees, formula=form, data=imports85, search="exhaustive", method="anova", split="MSE",
                                mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_imports85<-growRF_Parallel(ntrees=ntrees, formula=form, data=imports85, search="sss", method="anova", split="MSE",
                                 mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_imports85<-growRF_Parallel(ntrees=ntrees, formula=form, data=imports85, search="exhaustive", method="anova", split="MSE",
                                mtry=8, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_imports85<-growRF_Parallel(ntrees=ntrees, formula=form, data=imports85, search="ar", method="anova", split="MSE",
                                mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  mseRF_imports85[sim] <- mean((predictRF(rf_imports85,imports85,checkCases=TRUE)-imports85$Response)^2)
  mseSSS_imports85[sim] <- mean((predictRF(sss_imports85,imports85,checkCases=TRUE)-imports85$Response)^2)
  mseER_imports85[sim] <- mean((predictRF(er_imports85,imports85,checkCases=TRUE)-imports85$Response)^2)
  mseAR_imports85[sim] <- mean((predictRF(ar_imports85,imports85,checkCases=TRUE)-imports85$Response)^2)
}

outData <- data.frame(dataset=rep("Imports85", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_imports85, mseSSS_imports85, mseER_imports85, mseAR_imports85))

#write.table(outData, file = "Results/predImports85.csv", sep = ",", row.names=FALSE)





