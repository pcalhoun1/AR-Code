#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="contData.RData")

### Friedman1 ###

nsim=50

ntrees=100

# friedman1 #
friedman1<-contData$friedman1
form<-as.formula(paste0("Response ~ ",paste0("v",1:10,collapse=' + ')))

mseRF_friedman1 <- rep(NA, nsim); mseSSS_friedman1 <- rep(NA, nsim); mseER_friedman1 <- rep(NA, nsim); mseAR_friedman1 <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_friedman1<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman1, search="exhaustive", method="anova", split="MSE",
                                mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_friedman1<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman1, search="sss", method="anova", split="MSE",
                                 mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_friedman1<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman1, search="exhaustive", method="anova", split="MSE",
                                mtry=3, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_friedman1<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman1, search="ar", method="anova", split="MSE",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  mseRF_friedman1[sim] <- mean((predictRF(rf_friedman1,friedman1,checkCases=TRUE)-friedman1$Response)^2)
  mseSSS_friedman1[sim] <- mean((predictRF(sss_friedman1,friedman1,checkCases=TRUE)-friedman1$Response)^2)
  mseER_friedman1[sim] <- mean((predictRF(er_friedman1,friedman1,checkCases=TRUE)-friedman1$Response)^2)
  mseAR_friedman1[sim] <- mean((predictRF(ar_friedman1,friedman1,checkCases=TRUE)-friedman1$Response)^2)
}

outData <- data.frame(dataset=rep("Friedman1", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_friedman1, mseSSS_friedman1, mseER_friedman1, mseAR_friedman1))

write.table(outData, file = "predFriedman1_subsample.csv", sep = ",", row.names=FALSE)





