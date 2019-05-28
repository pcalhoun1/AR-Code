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

misclassAR_german <- rep(NA, 6*nsim)
aucAR_german <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_german<-growRF_Parallel(ntrees=ntrees, formula=form, data=german, search="ar", method="class", split="gini",
                               mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    
    lvls<-levels(german$Class)
    predAR <- predictRF(ar_german,german,checkCases=TRUE)
    misclassAR_german[index] <- mean(lvls[1+(predAR>0.5)] != german$Class)
    aucAR_german[index] <- auc(german$Class, predAR, direction="<")[1]
    index <- index + 1
  }
}

outData <- data.frame(dataset=rep("German", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      misclass = misclassAR_german, auc = aucAR_german)

#write.table(outData, file = "Results/minpvaluePredGerman.csv", sep = ",", row.names=FALSE)





