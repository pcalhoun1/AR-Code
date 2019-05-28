#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="binData.RData")

### Ringnorm ###

nsim=50

ntrees=1000


# ringnorm #
ringnorm<-binData$ringnorm
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:20),collapse=' + ')))

misclassRF_ringnorm <- rep(NA, nsim); misclassSSS_ringnorm <- rep(NA, nsim); misclassER_ringnorm <- rep(NA, nsim); misclassAR_ringnorm <- rep(NA, nsim)
aucRF_ringnorm <- rep(NA, nsim); aucSSS_ringnorm <- rep(NA, nsim); aucER_ringnorm <- rep(NA, nsim); aucAR_ringnorm <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_ringnorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=ringnorm, search="exhaustive", method="class", split="gini",
                               mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_ringnorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=ringnorm, search="sss", method="class", split="gini",
                                mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_ringnorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=ringnorm, search="exhaustive", method="class", split="gini",
                               mtry=4, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_ringnorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=ringnorm, search="ar", method="class", split="gini",
                               mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  lvls<-levels(ringnorm$Class)
  
  predRF <- predictRF(rf_ringnorm,ringnorm,checkCases=TRUE)
  predSSS <- predictRF(sss_ringnorm,ringnorm,checkCases=TRUE)
  predER <- predictRF(er_ringnorm,ringnorm,checkCases=TRUE)
  predAR <- predictRF(ar_ringnorm,ringnorm,checkCases=TRUE)
  
  misclassRF_ringnorm[sim] <- mean(lvls[1+(predRF>0.5)] != ringnorm$Class)
  misclassSSS_ringnorm[sim] <- mean(lvls[1+(predSSS>0.5)] != ringnorm$Class)
  misclassER_ringnorm[sim] <- mean(lvls[1+(predER>0.5)] != ringnorm$Class)
  misclassAR_ringnorm[sim] <- mean(lvls[1+(predAR>0.5)] != ringnorm$Class)
  
  aucRF_ringnorm[sim] <- auc(ringnorm$Class, predRF, direction="<")[1]
  aucSSS_ringnorm[sim] <- auc(ringnorm$Class, predSSS, direction="<")[1]
  aucER_ringnorm[sim] <- auc(ringnorm$Class, predER, direction="<")[1]
  aucAR_ringnorm[sim] <- auc(ringnorm$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Ringnorm", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_ringnorm, misclassSSS_ringnorm, misclassER_ringnorm, misclassAR_ringnorm),
                      auc = c(aucRF_ringnorm, aucSSS_ringnorm, aucER_ringnorm, aucAR_ringnorm))

write.table(outData, file = "predRingnorm_subsample.csv", sep = ",", row.names=FALSE)





