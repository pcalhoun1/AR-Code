#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('RF functions 20JAN19.R')
load(file="contData.RData")

### Abalone ###

nsim=50

ntrees=100

# abalone #
abalone<-contData$abalone
form<-as.formula("Response ~ sex + length + diameter + height + weight.w + weight.s + weight.v + weight.sh")

mseRF_abalone <- rep(NA, nsim); mseSSS_abalone <- rep(NA, nsim); mseER_abalone <- rep(NA, nsim); mseAR_abalone <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_abalone<-growRF_Parallel(ntrees=ntrees, formula=form, data=abalone, search="exhaustive", method="anova", split="MSE",
                              mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', useRpart=TRUE, iseed=sim)
  sss_abalone<-growRF_Parallel(ntrees=ntrees, formula=form, data=abalone, search="sss", method="anova", split="MSE",
                               mtry=3, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='subsample', iseed=sim)
  er_abalone<-growRF_Parallel(ntrees=ntrees, formula=form, data=abalone, search="exhaustive", method="anova", split="MSE",
                              mtry=3, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='subsample', iseed=sim)
  ar_abalone<-growRF_Parallel(ntrees=ntrees, formula=form, data=abalone, search="ar", method="anova", split="MSE",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='subsample', iseed=sim)
  
  mseRF_abalone[sim] <- mean((predictRF(rf_abalone,abalone,checkCases=TRUE)-abalone$Response)^2)
  mseSSS_abalone[sim] <- mean((predictRF(sss_abalone,abalone,checkCases=TRUE)-abalone$Response)^2)
  mseER_abalone[sim] <- mean((predictRF(er_abalone,abalone,checkCases=TRUE)-abalone$Response)^2)
  mseAR_abalone[sim] <- mean((predictRF(ar_abalone,abalone,checkCases=TRUE)-abalone$Response)^2)
}

outData <- data.frame(dataset=rep("Abalone", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_abalone, mseSSS_abalone, mseER_abalone, mseAR_abalone))

write.table(outData, file = "predAbalone_subsample.csv", sep = ",", row.names=FALSE)





