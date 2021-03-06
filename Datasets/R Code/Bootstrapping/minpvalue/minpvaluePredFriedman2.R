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

mseAR_friedman2 <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_friedman2<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman2, search="ar", method="anova", split="MSE",
                                  mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_friedman2[index] <- mean((predictRF(ar_friedman2,friedman2,checkCases=TRUE)-friedman2$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Friedman2", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_friedman2)

#write.table(outData, file = "Results/minpvaluePredFriedman2.csv", sep = ",", row.names=FALSE)





