#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")

### German ###

nsim=50

ntrees=1000


# german #
german<-binData$german
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:20),collapse=' + ')))

misclassRF_german <- rep(NA, nsim); misclassSSS_german <- rep(NA, nsim); misclassER_german <- rep(NA, nsim); misclassAR_german <- rep(NA, nsim)
aucRF_german <- rep(NA, nsim); aucSSS_german <- rep(NA, nsim); aucER_german <- rep(NA, nsim); aucAR_german <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_german<-growRF_Parallel(ntrees=ntrees, formula=form, data=german, search="exhaustive", method="class", split="gini",
                             mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_german<-growRF_Parallel(ntrees=ntrees, formula=form, data=german, search="sss", method="class", split="gini",
                              mtry=4, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_german<-growRF_Parallel(ntrees=ntrees, formula=form, data=german, search="exhaustive", method="class", split="gini",
                             mtry=4, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_german<-growRF_Parallel(ntrees=ntrees, formula=form, data=german, search="ar", method="class", split="gini",
                             mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  lvls<-levels(german$Class)
  
  predRF <- predictRF(rf_german,german,checkCases=TRUE)
  predSSS <- predictRF(sss_german,german,checkCases=TRUE)
  predER <- predictRF(er_german,german,checkCases=TRUE)
  predAR <- predictRF(ar_german,german,checkCases=TRUE)
  
  misclassRF_german[sim] <- mean(lvls[1+(predRF>0.5)] != german$Class)
  misclassSSS_german[sim] <- mean(lvls[1+(predSSS>0.5)] != german$Class)
  misclassER_german[sim] <- mean(lvls[1+(predER>0.5)] != german$Class)
  misclassAR_german[sim] <- mean(lvls[1+(predAR>0.5)] != german$Class)
  
  aucRF_german[sim] <- auc(german$Class, predRF, direction="<")[1]
  aucSSS_german[sim] <- auc(german$Class, predSSS, direction="<")[1]
  aucER_german[sim] <- auc(german$Class, predER, direction="<")[1]
  aucAR_german[sim] <- auc(german$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("German", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_german, misclassSSS_german, misclassER_german, misclassAR_german),
                      auc = c(aucRF_german, aucSSS_german, aucER_german, aucAR_german))

#write.table(outData, file = "Results/predGerman.csv", sep = ",", row.names=FALSE)





