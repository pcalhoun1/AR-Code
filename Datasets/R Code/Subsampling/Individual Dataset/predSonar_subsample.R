#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="binData.RData")

### Sonar ###

nsim=50

ntrees=1000


# sonar #
sonar<-binData$sonar
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:60),collapse=' + ')))

misclassRF_sonar <- rep(NA, nsim); misclassSSS_sonar <- rep(NA, nsim); misclassER_sonar <- rep(NA, nsim); misclassAR_sonar <- rep(NA, nsim)
aucRF_sonar <- rep(NA, nsim); aucSSS_sonar <- rep(NA, nsim); aucER_sonar <- rep(NA, nsim); aucAR_sonar <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_sonar<-growRF_Parallel(ntrees=ntrees, formula=form, data=sonar, search="exhaustive", method="class", split="gini",
                            mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_sonar<-growRF_Parallel(ntrees=ntrees, formula=form, data=sonar, search="sss", method="class", split="gini",
                             mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_sonar<-growRF_Parallel(ntrees=ntrees, formula=form, data=sonar, search="exhaustive", method="class", split="gini",
                            mtry=8, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_sonar<-growRF_Parallel(ntrees=ntrees, formula=form, data=sonar, search="ar", method="class", split="gini",
                            mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  lvls<-levels(sonar$Class)
  
  predRF <- predictRF(rf_sonar,sonar,checkCases=TRUE)
  predSSS <- predictRF(sss_sonar,sonar,checkCases=TRUE)
  predER <- predictRF(er_sonar,sonar,checkCases=TRUE)
  predAR <- predictRF(ar_sonar,sonar,checkCases=TRUE)
  
  misclassRF_sonar[sim] <- mean(lvls[1+(predRF>0.5)] != sonar$Class)
  misclassSSS_sonar[sim] <- mean(lvls[1+(predSSS>0.5)] != sonar$Class)
  misclassER_sonar[sim] <- mean(lvls[1+(predER>0.5)] != sonar$Class)
  misclassAR_sonar[sim] <- mean(lvls[1+(predAR>0.5)] != sonar$Class)
  
  aucRF_sonar[sim] <- auc(sonar$Class, predRF, direction="<")[1]
  aucSSS_sonar[sim] <- auc(sonar$Class, predSSS, direction="<")[1]
  aucER_sonar[sim] <- auc(sonar$Class, predER, direction="<")[1]
  aucAR_sonar[sim] <- auc(sonar$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Sonar", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_sonar, misclassSSS_sonar, misclassER_sonar, misclassAR_sonar),
                      auc = c(aucRF_sonar, aucSSS_sonar, aucER_sonar, aucAR_sonar))

write.table(outData, file = "predSonar_subsample.csv", sep = ",", row.names=FALSE)





