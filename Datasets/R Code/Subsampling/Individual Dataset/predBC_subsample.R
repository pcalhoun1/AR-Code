#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="binData.RData")

### Breast Cancer ###

nsim=50

ntrees=1000


# breastCancer #
breastCancer <- binData$breastCancer
form <- as.formula("Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion + Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses")

misclassRF_BC <- rep(NA, nsim); misclassSSS_BC <- rep(NA, nsim); misclassER_BC <- rep(NA, nsim); misclassAR_BC <- rep(NA, nsim)
aucRF_BC <- rep(NA, nsim); aucSSS_BC <- rep(NA, nsim); aucER_BC <- rep(NA, nsim); aucAR_BC <- rep(NA, nsim)

for (sim in 1:nsim) {

  rf_breastCancer <- growRF_Parallel(ntrees=ntrees, formula=form, data=breastCancer, search="exhaustive", method="class", split="gini",
                                     mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_breastCancer <- growRF_Parallel(ntrees=ntrees, formula=form, data=breastCancer, search="sss", method="class", split="gini",
                                      mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_breastCancer <- growRF_Parallel(ntrees=ntrees, formula=form, data=breastCancer, search="exhaustive", method="class", split="gini",
                                     mtry=3, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_breastCancer <- growRF_Parallel(ntrees=ntrees, formula=form, data=breastCancer, search="ar", method="class", split="gini",
                                     mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  lvls <- levels(breastCancer$Class)
  
  predRF <- predictRF(rf_breastCancer,breastCancer,checkCases=TRUE)
  predSSS <- predictRF(sss_breastCancer,breastCancer,checkCases=TRUE)
  predER <- predictRF(er_breastCancer,breastCancer,checkCases=TRUE)
  predAR <- predictRF(ar_breastCancer,breastCancer,checkCases=TRUE)
  
  misclassRF_BC[sim] <- mean(lvls[1+(predRF>0.5)] != breastCancer$Class)
  misclassSSS_BC[sim] <- mean(lvls[1+(predSSS>0.5)] != breastCancer$Class)
  misclassER_BC[sim] <- mean(lvls[1+(predER>0.5)] != breastCancer$Class)
  misclassAR_BC[sim] <- mean(lvls[1+(predAR>0.5)] != breastCancer$Class)
  
  aucRF_BC[sim] <- auc(breastCancer$Class, predRF, direction="<")[1]
  aucSSS_BC[sim] <- auc(breastCancer$Class, predSSS, direction="<")[1]
  aucER_BC[sim] <- auc(breastCancer$Class, predER, direction="<")[1]
  aucAR_BC[sim] <- auc(breastCancer$Class, predAR, direction="<")[1]
}

outData <- data.frame(dataset=rep("Breast Cancer", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      misclass = c(misclassRF_BC, misclassSSS_BC, misclassER_BC, misclassAR_BC),
                      auc = c(aucRF_BC, aucSSS_BC, aucER_BC, aucAR_BC))

write.table(outData, file = "predBC_subsample.csv", sep = ",", row.names=FALSE)





