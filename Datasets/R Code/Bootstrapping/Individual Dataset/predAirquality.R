#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")

### Airquality ###

nsim=50

ntrees=100

# airquality #
airquality<-contData$airquality
form<-as.formula("Response ~ Solar.R + Wind + Temp + Month + Day")

mseRF_airquality <- rep(NA, nsim); mseSSS_airquality <- rep(NA, nsim); mseER_airquality <- rep(NA, nsim); mseAR_airquality <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_airquality<-growRF_Parallel(ntrees=ntrees, formula=form, data=airquality, search="exhaustive", method="anova", split="MSE",
                                 mtry=2, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_airquality<-growRF_Parallel(ntrees=ntrees, formula=form, data=airquality, search="sss", method="anova", split="MSE",
                                  mtry=2, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_airquality<-growRF_Parallel(ntrees=ntrees, formula=form, data=airquality, search="exhaustive", method="anova", split="MSE",
                                 mtry=2, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_airquality<-growRF_Parallel(ntrees=ntrees, formula=form, data=airquality, search="ar", method="anova", split="MSE",
                                 mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  mseRF_airquality[sim] <- mean((predictRF(rf_airquality,airquality,checkCases=TRUE)-airquality$Response)^2)
  mseSSS_airquality[sim] <- mean((predictRF(sss_airquality,airquality,checkCases=TRUE)-airquality$Response)^2)
  mseER_airquality[sim] <- mean((predictRF(er_airquality,airquality,checkCases=TRUE)-airquality$Response)^2)
  mseAR_airquality[sim] <- mean((predictRF(ar_airquality,airquality,checkCases=TRUE)-airquality$Response)^2)
}

outData <- data.frame(dataset=rep("Airquality", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_airquality, mseSSS_airquality, mseER_airquality, mseAR_airquality))

#write.table(outData, file = "Results/predAirquality.csv", sep = ",", row.names=FALSE)





