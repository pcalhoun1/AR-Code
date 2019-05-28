#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="binData.RData")

### Twonorm ###

nsim=50

ntrees=1000


# twonorm #
twonorm<-binData$twonorm
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:20),collapse=' + ')))

misclassRF_twonorm <- rep(NA, nsim); misclassSSS_twonorm <- rep(NA, nsim); misclassER_twonorm <- rep(NA, nsim); misclassAR_twonorm <- rep(NA, nsim)
aucRF_twonorm <- rep(NA, nsim); aucSSS_twonorm <- rep(NA, nsim); aucER_twonorm <- rep(NA, nsim); aucAR_twonorm <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_twonorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=twonorm, search="exhaustive", method="class", split="gini",
                              mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_twonorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=twonorm, search="sss", method="class", split="gini",
                               mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_twonorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=twonorm, search="exhaustive", method="class", split="gini",
                              mtry=4, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_twonorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=twonorm, search="ar", method="class", split="gini",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  lvls<-levels(twonorm$Class)
  
  predRF <- predictRF(rf_twonorm,twonorm,checkCases=TRUE)
  predSSS <- predictRF(sss_twonorm,twonorm,checkCases=TRUE)
  predER <- predictRF(er_twonorm,twonorm,checkCases=TRUE)
  predAR <- predictRF(ar_twonorm,twonorm,checkCases=TRUE)
  
  misclassRF_twonorm[sim] <- mean(lvls[1+(predRF>0.5)] != twonorm$Class)
  misclassSSS_twonorm[sim] <- mean(lvls[1+(predSSS>0.5)] != twonorm$Class)
  misclassER_twonorm[sim] <- mean(lvls[1+(predER>0.5)] != twonorm$Class)
  misclassAR_twonorm[sim] <- mean(lvls[1+(predAR>0.5)] != twonorm$Class)
  
  aucRF_twonorm[sim] <- auc(twonorm$Class, predRF, direction="<")[1]
  aucSSS_twonorm[sim] <- auc(twonorm$Class, predSSS, direction="<")[1]
  aucER_twonorm[sim] <- auc(twonorm$Class, predER, direction="<")[1]
  aucAR_twonorm[sim] <- auc(twonorm$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Twonorm", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_twonorm, misclassSSS_twonorm, misclassER_twonorm, misclassAR_twonorm),
                      auc = c(aucRF_twonorm, aucSSS_twonorm, aucER_twonorm, aucAR_twonorm))

write.table(outData, file = "predTwonorm_subsample.csv", sep = ",", row.names=FALSE)





