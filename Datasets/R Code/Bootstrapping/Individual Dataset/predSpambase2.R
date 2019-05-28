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

misclassER_spambase <- rep(NA, nsim); misclassAR_spambase <- rep(NA, nsim)
aucER_spambase <- rep(NA, nsim); aucAR_spambase <- rep(NA, nsim)

for (sim in 1:nsim) {
  er_spambase<-growRF_Parallel(ntrees=ntrees, formula=form, data=spambase, search="exhaustive", method="class", split="gini",
                               mtry=8, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_spambase<-growRF_Parallel(ntrees=ntrees, formula=form, data=spambase, search="ar", method="class", split="gini",
                               mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  lvls<-levels(spambase$Class)
  
  predER <- predictRF(er_spambase,spambase,checkCases=TRUE)
  predAR <- predictRF(ar_spambase,spambase,checkCases=TRUE)
  
  misclassER_spambase[sim] <- mean(lvls[1+(predER>0.5)] != spambase$Class)
  misclassAR_spambase[sim] <- mean(lvls[1+(predAR>0.5)] != spambase$Class)
  
  aucER_spambase[sim] <- auc(spambase$Class, predER, direction="<")[1]
  aucAR_spambase[sim] <- auc(spambase$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Spambase", 2*nsim), sim=rep(1:nsim, 2), method=rep(c("ER", "AR"), each=nsim),
                      misclass = c(misclassER_spambase, misclassAR_spambase),
                      auc = c(aucER_spambase, aucAR_spambase))

#write.table(outData, file = "Results/predSpambase2.csv", sep = ",", row.names=FALSE)





