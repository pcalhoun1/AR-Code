#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")


### Friedman1 ###

nsim=50

ntrees=100

# friedman1 #
friedman1<-contData$friedman1
form<-as.formula(paste0("Response ~ ",paste0("v",1:10,collapse=' + ')))

mseAR_friedman1 <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_friedman1<-growRF_Parallel(ntrees=ntrees, formula=form, data=friedman1, search="ar", method="anova", split="MSE",
                                  mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_friedman1[index] <- mean((predictRF(ar_friedman1,friedman1,checkCases=TRUE)-friedman1$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Friedman1", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_friedman1)

#write.table(outData, file = "Results/minpvaluePredFriedman1.csv", sep = ",", row.names=FALSE)





