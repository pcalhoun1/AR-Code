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

misclassRF_spambase <- rep(NA, nsim); misclassSSS_spambase <- rep(NA, nsim)
aucRF_spambase <- rep(NA, nsim); aucSSS_spambase <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_spambase<-growRF_Parallel(ntrees=ntrees, formula=form, data=spambase, search="exhaustive", method="class", split="gini",
                               mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_spambase<-growRF_Parallel(ntrees=ntrees, formula=form, data=spambase, search="sss", method="class", split="gini",
                                mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)

  lvls<-levels(spambase$Class)
  
  predRF <- predictRF(rf_spambase,spambase,checkCases=TRUE)
  predSSS <- predictRF(sss_spambase,spambase,checkCases=TRUE)

  misclassRF_spambase[sim] <- mean(lvls[1+(predRF>0.5)] != spambase$Class)
  misclassSSS_spambase[sim] <- mean(lvls[1+(predSSS>0.5)] != spambase$Class)

  aucRF_spambase[sim] <- auc(spambase$Class, predRF, direction="<")[1]
  aucSSS_spambase[sim] <- auc(spambase$Class, predSSS, direction="<")[1]
}

outData <- data.frame(dataset=rep("Spambase", 2*nsim), sim=rep(1:nsim, 2), method=rep(c("RF", "SSS"), each=nsim),
                      misclass = c(misclassRF_spambase, misclassSSS_spambase),
                      auc = c(aucRF_spambase, aucSSS_spambase))

#write.table(outData, file = "Results/predSpambase1.csv", sep = ",", row.names=FALSE)





