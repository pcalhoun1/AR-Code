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

mseAR_airquality <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_airquality<-growRF_Parallel(ntrees=ntrees, formula=form, data=airquality, search="ar", method="anova", split="MSE",
                                   mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_airquality[index] <- mean((predictRF(ar_airquality,airquality,checkCases=TRUE)-airquality$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Airquality", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_airquality)

#write.table(outData, file = "Results/minpvaluePredAirquality.csv", sep = ",", row.names=FALSE)





