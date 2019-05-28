#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="contData.RData")

### Ailerons ###

nsim=50

ntrees=100

# ailerons #
ailerons<-contData$ailerons
form<-as.formula(paste0("Response ~ ",paste0("V",1:40,collapse=' + ')))

mseRF_ailerons <- rep(NA, nsim); mseSSS_ailerons <- rep(NA, nsim); mseER_ailerons <- rep(NA, nsim); mseAR_ailerons <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_ailerons<-growRF_Parallel(ntrees=ntrees, formula=form, data=ailerons, search="exhaustive", method="anova", split="MSE",
                               mtry=13, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_ailerons<-growRF_Parallel(ntrees=ntrees, formula=form, data=ailerons, search="sss", method="anova", split="MSE",
                                mtry=13, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_ailerons<-growRF_Parallel(ntrees=ntrees, formula=form, data=ailerons, search="exhaustive", method="anova", split="MSE",
                               mtry=13, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_ailerons<-growRF_Parallel(ntrees=ntrees, formula=form, data=ailerons, search="ar", method="anova", split="MSE",
                               mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  mseRF_ailerons[sim] <- mean((predictRF(rf_ailerons,ailerons,checkCases=TRUE)-ailerons$Response)^2)
  mseSSS_ailerons[sim] <- mean((predictRF(sss_ailerons,ailerons,checkCases=TRUE)-ailerons$Response)^2)
  mseER_ailerons[sim] <- mean((predictRF(er_ailerons,ailerons,checkCases=TRUE)-ailerons$Response)^2)
  mseAR_ailerons[sim] <- mean((predictRF(ar_ailerons,ailerons,checkCases=TRUE)-ailerons$Response)^2)
}

outData <- data.frame(dataset=rep("Ailerons", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_ailerons, mseSSS_ailerons, mseER_ailerons, mseAR_ailerons))

write.table(outData, file = "predAilerons_subsample.csv", sep = ",", row.names=FALSE)





