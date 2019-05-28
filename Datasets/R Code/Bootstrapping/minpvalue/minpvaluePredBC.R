#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")


### Breast Cancer ###

nsim=50

ntrees=1000


# breastCancer #
breastCancer <- binData$breastCancer
form <- as.formula("Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion + Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses")

misclassAR_BC <- rep(NA, 6*nsim)
aucAR_BC <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_breastCancer <- growRF_Parallel(ntrees=ntrees, formula=form, data=breastCancer, search="ar", method="class", split="gini",
                                       mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls <- levels(breastCancer$Class)
    predAR <- predictRF(ar_breastCancer,breastCancer,checkCases=TRUE)
    misclassAR_BC[index] <- mean(lvls[1+(predAR>0.5)] != breastCancer$Class)
    aucAR_BC[index] <- auc(breastCancer$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Breast Cancer", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_BC, auc = aucAR_BC)

#write.table(outData, file = "Results/minpvaluePredBC.csv", sep = ",", row.names=FALSE)





