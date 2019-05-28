#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")

### Friedman3 ###

nsim=50

ntrees=100

# friedman3 #
friedman3<-contData$friedman3
form<-as.formula(paste0("Response ~ ",paste0("v",1:4,collapse=' + ')))

mseRF_friedman3 <- rep(NA, nsim); mseSSS_friedman3 <- rep(NA, nsim); mseER_friedman3 <- rep(NA, nsim); mseAR_friedman3 <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_friedman3<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman3, search="exhaustive", method="anova", split="MSE",
                                mtry=1, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_friedman3<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman3, search="sss", method="anova", split="MSE",
                                 mtry=1, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_friedman3<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman3, search="exhaustive", method="anova", split="MSE",
                                mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_friedman3<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman3, search="ar", method="anova", split="MSE",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  mseRF_friedman3[sim] <- mean((predictRF(rf_friedman3,friedman3,checkCases=TRUE)-friedman3$Response)^2)
  mseSSS_friedman3[sim] <- mean((predictRF(sss_friedman3,friedman3,checkCases=TRUE)-friedman3$Response)^2)
  mseER_friedman3[sim] <- mean((predictRF(er_friedman3,friedman3,checkCases=TRUE)-friedman3$Response)^2)
  mseAR_friedman3[sim] <- mean((predictRF(ar_friedman3,friedman3,checkCases=TRUE)-friedman3$Response)^2)
}

outData <- data.frame(dataset=rep("Friedman3", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_friedman3, mseSSS_friedman3, mseER_friedman3, mseAR_friedman3))

#write.table(outData, file = "Results/predFriedman3.csv", sep = ",", row.names=FALSE)





