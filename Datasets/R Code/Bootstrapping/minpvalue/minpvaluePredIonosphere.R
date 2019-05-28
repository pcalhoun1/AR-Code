#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")


### Ionosphere ###

nsim=50

ntrees=1000

# ionosphere #
ionosphere<-binData$ionosphere
form<-as.formula(paste0("Class ~ V1 + ",paste(paste0("V",3:34),collapse=' + ')))

misclassAR_ionosphere <- rep(NA, 6*nsim)
aucAR_ionosphere <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_ionosphere<-growRF_Parallel(ntrees=ntrees, formula=form, data=ionosphere, search="ar", method="class", split="gini",
                                   mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(ionosphere$Class)
    predAR <- predictRF(ar_ionosphere,ionosphere,checkCases=TRUE)
    misclassAR_ionosphere[index] <- mean(lvls[1+(predAR>0.5)] != ionosphere$Class)
    aucAR_ionosphere[index] <- auc(ionosphere$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Ionosphere", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_ionosphere, auc = aucAR_ionosphere)

#write.table(outData, file = "Results/minpvaluePredIonosphere.csv", sep = ",", row.names=FALSE)





