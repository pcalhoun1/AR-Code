#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")

### Birthwt ###

nsim=50

ntrees=1000


# birthwt #
birthwt<-binData$birthwt
form<-as.formula("Class ~ age + lwt + race + smoke + ptl + ht + ui + ftv")

misclassRF_birthwt <- rep(NA, nsim); misclassSSS_birthwt <- rep(NA, nsim); misclassER_birthwt <- rep(NA, nsim); misclassAR_birthwt <- rep(NA, nsim)
aucRF_birthwt <- rep(NA, nsim); aucSSS_birthwt <- rep(NA, nsim); aucER_birthwt <- rep(NA, nsim); aucAR_birthwt <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_birthwt<-growRF_Parallel(ntrees=ntrees, formula=form, data=birthwt, search="exhaustive", method="class", split="gini",
                              mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_birthwt<-growRF_Parallel(ntrees=ntrees, formula=form, data=birthwt, search="sss", method="class", split="gini",
                               mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_birthwt<-growRF_Parallel(ntrees=ntrees, formula=form, data=birthwt, search="exhaustive", method="class", split="gini",
                              mtry=3, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_birthwt<-growRF_Parallel(ntrees=ntrees, formula=form, data=birthwt, search="ar", method="class", split="gini",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  lvls<-levels(birthwt$Class)
  
  predRF <- predictRF(rf_birthwt,birthwt,checkCases=TRUE)
  predSSS <- predictRF(sss_birthwt,birthwt,checkCases=TRUE)
  predER <- predictRF(er_birthwt,birthwt,checkCases=TRUE)
  predAR <- predictRF(ar_birthwt,birthwt,checkCases=TRUE)
  
  misclassRF_birthwt[sim] <- mean(lvls[1+(predRF>0.5)] != birthwt$Class)
  misclassSSS_birthwt[sim] <- mean(lvls[1+(predSSS>0.5)] != birthwt$Class)
  misclassER_birthwt[sim] <- mean(lvls[1+(predER>0.5)] != birthwt$Class)
  misclassAR_birthwt[sim] <- mean(lvls[1+(predAR>0.5)] != birthwt$Class)
  
  aucRF_birthwt[sim] <- auc(birthwt$Class, predRF, direction="<")[1]
  aucSSS_birthwt[sim] <- auc(birthwt$Class, predSSS, direction="<")[1]
  aucER_birthwt[sim] <- auc(birthwt$Class, predER, direction="<")[1]
  aucAR_birthwt[sim] <- auc(birthwt$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Birthwt", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_birthwt, misclassSSS_birthwt, misclassER_birthwt, misclassAR_birthwt),
                      auc = c(aucRF_birthwt, aucSSS_birthwt, aucER_birthwt, aucAR_birthwt))

#write.table(outData, file = "Results/predBirthwt.csv", sep = ",", row.names=FALSE)





