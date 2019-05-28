#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../Data/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")


### Twonorm ###

nsim=50

ntrees=1000

# twonorm #
twonorm<-binData$twonorm
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:20),collapse=' + ')))

misclassAR_twonorm <- rep(NA, 6*nsim)
aucAR_twonorm <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_twonorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=twonorm, search="ar", method="class", split="gini",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(twonorm$Class)
    predAR <- predictRF(ar_twonorm,twonorm,checkCases=TRUE)
    misclassAR_twonorm[index] <- mean(lvls[1+(predAR>0.5)] != twonorm$Class)
    aucAR_twonorm[index] <- auc(twonorm$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Twonorm", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_twonorm, auc = aucAR_twonorm)

#write.table(outData, file = "Results/minpvaluePredTwonorm.csv", sep = ",", row.names=FALSE)





