#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")

### Liver ###

nsim=50

ntrees=1000


# liver #
liver<-binData$liver
form<-as.formula("Class ~ mcv + alkphos + sgpt + sgot + gammagt + drinks")

misclassRF_liver <- rep(NA, nsim); misclassSSS_liver <- rep(NA, nsim); misclassER_liver <- rep(NA, nsim); misclassAR_liver <- rep(NA, nsim)
aucRF_liver <- rep(NA, nsim); aucSSS_liver <- rep(NA, nsim); aucER_liver <- rep(NA, nsim); aucAR_liver <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_liver<-growRF_Parallel(ntrees=ntrees, formula=form, data=liver, search="exhaustive", method="class", split="gini",
                            mtry=2, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_liver<-growRF_Parallel(ntrees=ntrees, formula=form, data=liver, search="sss", method="class", split="gini",
                             mtry=2, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_liver<-growRF_Parallel(ntrees=ntrees, formula=form, data=liver, search="exhaustive", method="class", split="gini",
                            mtry=2, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_liver<-growRF_Parallel(ntrees=ntrees, formula=form, data=liver, search="ar", method="class", split="gini",
                            mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  lvls<-levels(liver$Class)
  
  predRF <- predictRF(rf_liver,liver,checkCases=TRUE)
  predSSS <- predictRF(sss_liver,liver,checkCases=TRUE)
  predER <- predictRF(er_liver,liver,checkCases=TRUE)
  predAR <- predictRF(ar_liver,liver,checkCases=TRUE)
  
  misclassRF_liver[sim] <- mean(lvls[1+(predRF>0.5)] != liver$Class)
  misclassSSS_liver[sim] <- mean(lvls[1+(predSSS>0.5)] != liver$Class)
  misclassER_liver[sim] <- mean(lvls[1+(predER>0.5)] != liver$Class)
  misclassAR_liver[sim] <- mean(lvls[1+(predAR>0.5)] != liver$Class)
  
  aucRF_liver[sim] <- auc(liver$Class, predRF, direction="<")[1]
  aucSSS_liver[sim] <- auc(liver$Class, predSSS, direction="<")[1]
  aucER_liver[sim] <- auc(liver$Class, predER, direction="<")[1]
  aucAR_liver[sim] <- auc(liver$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Liver", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_liver, misclassSSS_liver, misclassER_liver, misclassAR_liver),
                      auc = c(aucRF_liver, aucSSS_liver, aucER_liver, aucAR_liver))

#write.table(outData, file = "Results/predLiver.csv", sep = ",", row.names=FALSE)





