#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")


### Elevators ###

nsim=50

ntrees=100

# elevators #
elevators<-contData$elevators
form<-as.formula(paste0("Response ~ ",paste0("V",1:18,collapse=' + ')))

mseAR_elevators <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_elevators<-growRF_Parallel(ntrees=ntrees, formula=form, data=elevators, search="ar", method="anova", split="MSE",
                                  mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_elevators[index] <- mean((predictRF(ar_elevators,elevators,checkCases=TRUE)-elevators$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Elevators", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_elevators)

#write.table(outData, file = "Results/minpvaluePredElevators.csv", sep = ",", row.names=FALSE)





