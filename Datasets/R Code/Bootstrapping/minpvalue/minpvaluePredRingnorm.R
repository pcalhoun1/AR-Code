#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/binData.RData")


### Ringnorm ###

nsim=50

ntrees=1000

# ringnorm #
ringnorm<-binData$ringnorm
form<-as.formula(paste0("Class ~ ",paste(paste0("V",1:20),collapse=' + ')))

misclassAR_ringnorm <- rep(NA, 6*nsim)
aucAR_ringnorm <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_ringnorm<-growRF_Parallel(ntrees=ntrees, formula=form, data=ringnorm, search="ar", method="class", split="gini",
                                 mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(ringnorm$Class)
    predAR <- predictRF(ar_ringnorm,ringnorm,checkCases=TRUE)
    misclassAR_ringnorm[index] <- mean(lvls[1+(predAR>0.5)] != ringnorm$Class)
    aucAR_ringnorm[index] <- auc(ringnorm$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("Ringnorm", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_ringnorm, auc = aucAR_ringnorm)

#write.table(outData, file = "Results/minpvaluePredRingnorm.csv", sep = ",", row.names=FALSE)





