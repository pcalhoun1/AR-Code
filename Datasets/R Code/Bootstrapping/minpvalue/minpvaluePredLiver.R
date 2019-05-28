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

misclassAR_liver <- rep(NA, 6*nsim)
aucAR_liver <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_liver<-growRF_Parallel(ntrees=ntrees, formula=form, data=liver, search="ar", method="class", split="gini",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(liver$Class)
    predAR <- predictRF(ar_liver,liver,checkCases=TRUE)
    misclassAR_liver[index] <- mean(lvls[1+(predAR>0.5)] != liver$Class)
    aucAR_liver[index] <- auc(liver$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Liver", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_liver, auc = aucAR_liver)

#write.table(outData, file = "Results/minpvaluePredLiver.csv", sep = ",", row.names=FALSE)





