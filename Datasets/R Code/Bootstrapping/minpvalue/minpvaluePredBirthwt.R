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

misclassAR_birthwt <- rep(NA, 6*nsim)
aucAR_birthwt <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_birthwt<-growRF_Parallel(ntrees=ntrees, formula=form, data=birthwt, search="ar", method="class", split="gini",
                                mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(birthwt$Class)
    predAR <- predictRF(ar_birthwt,birthwt,checkCases=TRUE)
    misclassAR_birthwt[index] <- mean(lvls[1+(predAR>0.5)] != birthwt$Class)
    aucAR_birthwt[index] <- auc(birthwt$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Birthwt", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_birthwt, auc = aucAR_birthwt)

#write.table(outData, file = "Results/minpvaluePredBirthwt.csv", sep = ",", row.names=FALSE)





