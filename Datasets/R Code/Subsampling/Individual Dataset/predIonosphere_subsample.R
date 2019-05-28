#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="binData.RData")

### Ionosphere ###

nsim=50

ntrees=1000


# ionosphere #
ionosphere<-binData$ionosphere
form<-as.formula(paste0("Class ~ V1 + ",paste(paste0("V",3:34),collapse=' + ')))

misclassRF_ionosphere <- rep(NA, nsim); misclassSSS_ionosphere <- rep(NA, nsim); misclassER_ionosphere <- rep(NA, nsim); misclassAR_ionosphere <- rep(NA, nsim)
aucRF_ionosphere <- rep(NA, nsim); aucSSS_ionosphere <- rep(NA, nsim); aucER_ionosphere <- rep(NA, nsim); aucAR_ionosphere <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_ionosphere<-growRF_Parallel(ntrees=ntrees, formula=form, data=ionosphere, search="exhaustive", method="class", split="gini",
                                 mtry=6, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_ionosphere<-growRF_Parallel(ntrees=ntrees, formula=form, data=ionosphere, search="sss", method="class", split="gini",
                                  mtry=6, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_ionosphere<-growRF_Parallel(ntrees=ntrees, formula=form, data=ionosphere, search="exhaustive", method="class", split="gini",
                                 mtry=6, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_ionosphere<-growRF_Parallel(ntrees=ntrees, formula=form, data=ionosphere, search="ar", method="class", split="gini",
                                 mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  lvls<-levels(ionosphere$Class)
  
  predRF <- predictRF(rf_ionosphere,ionosphere,checkCases=TRUE)
  predSSS <- predictRF(sss_ionosphere,ionosphere,checkCases=TRUE)
  predER <- predictRF(er_ionosphere,ionosphere,checkCases=TRUE)
  predAR <- predictRF(ar_ionosphere,ionosphere,checkCases=TRUE)
  
  misclassRF_ionosphere[sim] <- mean(lvls[1+(predRF>0.5)] != ionosphere$Class)
  misclassSSS_ionosphere[sim] <- mean(lvls[1+(predSSS>0.5)] != ionosphere$Class)
  misclassER_ionosphere[sim] <- mean(lvls[1+(predER>0.5)] != ionosphere$Class)
  misclassAR_ionosphere[sim] <- mean(lvls[1+(predAR>0.5)] != ionosphere$Class)
  
  aucRF_ionosphere[sim] <- auc(ionosphere$Class, predRF, direction="<")[1]
  aucSSS_ionosphere[sim] <- auc(ionosphere$Class, predSSS, direction="<")[1]
  aucER_ionosphere[sim] <- auc(ionosphere$Class, predER, direction="<")[1]
  aucAR_ionosphere[sim] <- auc(ionosphere$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Ionosphere", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_ionosphere, misclassSSS_ionosphere, misclassER_ionosphere, misclassAR_ionosphere),
                      auc = c(aucRF_ionosphere, aucSSS_ionosphere, aucER_ionosphere, aucAR_ionosphere))

write.table(outData, file = "predIonosphere_subsample.csv", sep = ",", row.names=FALSE)





