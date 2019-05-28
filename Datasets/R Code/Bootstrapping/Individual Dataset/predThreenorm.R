#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")

### Threenorm ###

nsim=50

ntrees=1000


# threenorm #
threenorm<-binData$threenorm
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:20),collapse=' + ')))

misclassRF_threenorm <- rep(NA, nsim); misclassSSS_threenorm <- rep(NA, nsim); misclassER_threenorm <- rep(NA, nsim); misclassAR_threenorm <- rep(NA, nsim)
aucRF_threenorm <- rep(NA, nsim); aucSSS_threenorm <- rep(NA, nsim); aucER_threenorm <- rep(NA, nsim); aucAR_threenorm <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_threenorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=threenorm, search="exhaustive", method="class", split="gini",
                                mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_threenorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=threenorm, search="sss", method="class", split="gini",
                                 mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_threenorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=threenorm, search="exhaustive", method="class", split="gini",
                                mtry=4, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_threenorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=threenorm, search="ar", method="class", split="gini",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  lvls<-levels(threenorm$Class)
  
  predRF <- predictRF(rf_threenorm,threenorm,checkCases=TRUE)
  predSSS <- predictRF(sss_threenorm,threenorm,checkCases=TRUE)
  predER <- predictRF(er_threenorm,threenorm,checkCases=TRUE)
  predAR <- predictRF(ar_threenorm,threenorm,checkCases=TRUE)
  
  misclassRF_threenorm[sim] <- mean(lvls[1+(predRF>0.5)] != threenorm$Class)
  misclassSSS_threenorm[sim] <- mean(lvls[1+(predSSS>0.5)] != threenorm$Class)
  misclassER_threenorm[sim] <- mean(lvls[1+(predER>0.5)] != threenorm$Class)
  misclassAR_threenorm[sim] <- mean(lvls[1+(predAR>0.5)] != threenorm$Class)
  
  aucRF_threenorm[sim] <- auc(threenorm$Class, predRF, direction="<")[1]
  aucSSS_threenorm[sim] <- auc(threenorm$Class, predSSS, direction="<")[1]
  aucER_threenorm[sim] <- auc(threenorm$Class, predER, direction="<")[1]
  aucAR_threenorm[sim] <- auc(threenorm$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Threenorm", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_threenorm, misclassSSS_threenorm, misclassER_threenorm, misclassAR_threenorm),
                      auc = c(aucRF_threenorm, aucSSS_threenorm, aucER_threenorm, aucAR_threenorm))

#write.table(outData, file = "Results/predThreenorm.csv", sep = ",", row.names=FALSE)





