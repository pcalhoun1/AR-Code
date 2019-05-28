#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")


### Sonar ###

nsim=50

ntrees=1000

# sonar #
sonar<-binData$sonar
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:60),collapse=' + ')))

misclassAR_sonar <- rep(NA, 6*nsim)
aucAR_sonar <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_sonar<-growRF_Parallel(ntrees=ntrees, formula=form, data=sonar, search="ar", method="class", split="gini",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(sonar$Class)
    predAR <- predictRF(ar_sonar,sonar,checkCases=TRUE)
    misclassAR_sonar[index] <- mean(lvls[1+(predAR>0.5)] != sonar$Class)
    aucAR_sonar[index] <- auc(sonar$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Sonar", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_sonar, auc = aucAR_sonar)

#write.table(outData, file = "Results/minpvaluePredSonar.csv", sep = ",", row.names=FALSE)





