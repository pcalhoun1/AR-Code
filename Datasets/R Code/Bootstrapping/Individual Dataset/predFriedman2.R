#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")

### Friedman2 ###

nsim=50

ntrees=100

# friedman2 #
friedman2<-contData$friedman2
form<-as.formula(paste0("Response ~ ",paste0("v",1:4,collapse=' + ')))

mseRF_friedman2 <- rep(NA, nsim); mseSSS_friedman2 <- rep(NA, nsim); mseER_friedman2 <- rep(NA, nsim); mseAR_friedman2 <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_friedman2<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman2, search="exhaustive", method="anova", split="MSE",
                                mtry=1, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_friedman2<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman2, search="sss", method="anova", split="MSE",
                                 mtry=1, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_friedman2<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman2, search="exhaustive", method="anova", split="MSE",
                                mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_friedman2<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman2, search="ar", method="anova", split="MSE",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  mseRF_friedman2[sim] <- mean((predictRF(rf_friedman2,friedman2,checkCases=TRUE)-friedman2$Response)^2)
  mseSSS_friedman2[sim] <- mean((predictRF(sss_friedman2,friedman2,checkCases=TRUE)-friedman2$Response)^2)
  mseER_friedman2[sim] <- mean((predictRF(er_friedman2,friedman2,checkCases=TRUE)-friedman2$Response)^2)
  mseAR_friedman2[sim] <- mean((predictRF(ar_friedman2,friedman2,checkCases=TRUE)-friedman2$Response)^2)
}

outData <- data.frame(dataset=rep("Friedman2", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_friedman2, mseSSS_friedman2, mseER_friedman2, mseAR_friedman2))

#write.table(outData, file = "Results/predFriedman2.csv", sep = ",", row.names=FALSE)





