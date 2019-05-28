#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")


### Housing ###

nsim=50

ntrees=100

# housing #
housing<-contData$housing
form<-as.formula("Response ~ town + tract + lon + lat + crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + b + lstat")

mseAR_housing <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_housing<-growRF_Parallel(ntrees=ntrees, formula=form, data=housing, search="ar", method="anova", split="MSE",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_housing[index] <- mean((predictRF(ar_housing,housing,checkCases=TRUE)-housing$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Housing", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_housing)

#write.table(outData, file = "Results/minpvaluePredHousing.csv", sep = ",", row.names=FALSE)





