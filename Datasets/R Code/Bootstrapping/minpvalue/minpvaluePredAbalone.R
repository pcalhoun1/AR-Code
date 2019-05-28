#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")

### Abalone ###

nsim=50

ntrees=100

# abalone #
abalone<-contData$abalone
form<-as.formula("Response ~ sex + length + diameter + height + weight.w + weight.s + weight.v + weight.sh")

mseAR_abalone <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_abalone<-growRF_Parallel(ntrees=ntrees, formula=form, data=abalone, search="ar", method="anova", split="MSE",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_abalone[index] <- mean((predictRF(ar_abalone,abalone,checkCases=TRUE)-abalone$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Abalone", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_abalone)

#write.table(outData, file = "Results/minpvaluePredAbalone.csv", sep = ",", row.names=FALSE)





