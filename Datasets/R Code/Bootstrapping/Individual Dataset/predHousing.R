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

mseRF_housing <- rep(NA, nsim); mseSSS_housing <- rep(NA, nsim); mseER_housing <- rep(NA, nsim); mseAR_housing <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_housing<-growRF_Parallel(ntrees=ntrees, formula=form, data=housing, search="exhaustive", method="anova", split="MSE",
                              mtry=6, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_housing<-growRF_Parallel(ntrees=ntrees, formula=form, data=housing, search="sss", method="anova", split="MSE",
                               mtry=6, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_housing<-growRF_Parallel(ntrees=ntrees, formula=form, data=housing, search="exhaustive", method="anova", split="MSE",
                              mtry=6, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_housing<-growRF_Parallel(ntrees=ntrees, formula=form, data=housing, search="ar", method="anova", split="MSE",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  mseRF_housing[sim] <- mean((predictRF(rf_housing,housing,checkCases=TRUE)-housing$Response)^2)
  mseSSS_housing[sim] <- mean((predictRF(sss_housing,housing,checkCases=TRUE)-housing$Response)^2)
  mseER_housing[sim] <- mean((predictRF(er_housing,housing,checkCases=TRUE)-housing$Response)^2)
  mseAR_housing[sim] <- mean((predictRF(ar_housing,housing,checkCases=TRUE)-housing$Response)^2)
}

outData <- data.frame(dataset=rep("Housing", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_housing, mseSSS_housing, mseER_housing, mseAR_housing))

#write.table(outData, file = "Results/predHousing.csv", sep = ",", row.names=FALSE)





