#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="contData.RData")

### Elevators ###

nsim=50

ntrees=100

# elevators #
elevators<-contData$elevators
form<-as.formula(paste0("Response ~ ",paste0("V",1:18,collapse=' + ')))

mseRF_elevators <- rep(NA, nsim); mseSSS_elevators <- rep(NA, nsim); mseER_elevators <- rep(NA, nsim); mseAR_elevators <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_elevators<-growRF_Parallel(ntrees=ntrees, formula=form, data=elevators, search="exhaustive", method="anova", split="MSE",
                                mtry=6, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_elevators<-growRF_Parallel(ntrees=ntrees, formula=form, data=elevators, search="sss", method="anova", split="MSE",
                                 mtry=6, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_elevators<-growRF_Parallel(ntrees=ntrees, formula=form, data=elevators, search="exhaustive", method="anova", split="MSE",
                                mtry=6, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_elevators<-growRF_Parallel(ntrees=ntrees, formula=form, data=elevators, search="ar", method="anova", split="MSE",
                                mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  mseRF_elevators[sim] <- mean((predictRF(rf_elevators,elevators,checkCases=TRUE)-elevators$Response)^2)
  mseSSS_elevators[sim] <- mean((predictRF(sss_elevators,elevators,checkCases=TRUE)-elevators$Response)^2)
  mseER_elevators[sim] <- mean((predictRF(er_elevators,elevators,checkCases=TRUE)-elevators$Response)^2)
  mseAR_elevators[sim] <- mean((predictRF(ar_elevators,elevators,checkCases=TRUE)-elevators$Response)^2)
}

outData <- data.frame(dataset=rep("Elevators", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_elevators, mseSSS_elevators, mseER_elevators, mseAR_elevators))

write.table(outData, file = "predElevators_subsample.csv", sep = ",", row.names=FALSE)





