#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")


### Spambase ###

nsim=50

ntrees=1000

# spambase #
spambase<-binData$spambase
form<-as.formula(paste0("Class ~ ",paste0("V",1:57,collapse=" + ")))

misclassAR_spambase <- rep(NA, 6*nsim)
aucAR_spambase <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_spambase<-growRF_Parallel(ntrees=ntrees, formula=form, data=spambase, search="ar", method="class", split="gini",
                                 mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(spambase$Class)
    predAR <- predictRF(ar_spambase,spambase,checkCases=TRUE)
    misclassAR_spambase[index] <- mean(lvls[1+(predAR>0.5)] != spambase$Class)
    aucAR_spambase[index] <- auc(spambase$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Spambase", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_spambase, auc = aucAR_spambase)

#write.table(outData, file = "Results/minpvaluePredSpambase.csv", sep = ",", row.names=FALSE)





