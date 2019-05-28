#rm(list=ls(all=TRUE))

library(parallel)

getwd()
dir()

source('../../../../R Functions/RF functions 20JAN19.R')
load(file="../../../Data/contData.RData")

### Servo ###

nsim=50

ntrees=100

# servo #
servo<-contData$servo
form<-as.formula("Response ~ Motor + Screw + Pgain + Vgain")

mseRF_servo <- rep(NA, nsim); mseSSS_servo <- rep(NA, nsim); mseER_servo <- rep(NA, nsim); mseAR_servo <- rep(NA, nsim)

for (sim in 1:nsim) {
  rf_servo<-growRF_Parallel(ntrees=ntrees, formula=form, data=servo, search="exhaustive", method="anova", split="MSE",
                            mtry=1, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE, iseed=sim)
  sss_servo<-growRF_Parallel(ntrees=ntrees, formula=form, data=servo, search="sss", method="anova", split="MSE",
                             mtry=1, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap', iseed=sim)
  er_servo<-growRF_Parallel(ntrees=ntrees, formula=form, data=servo, search="exhaustive", method="anova", split="MSE",
                            mtry=1, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', iseed=sim)
  ar_servo<-growRF_Parallel(ntrees=ntrees, formula=form, data=servo, search="ar", method="anova", split="MSE",
                            mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap', iseed=sim)
  
  mseRF_servo[sim] <- mean((predictRF(rf_servo,servo,checkCases=TRUE)-servo$Response)^2)
  mseSSS_servo[sim] <- mean((predictRF(sss_servo,servo,checkCases=TRUE)-servo$Response)^2)
  mseER_servo[sim] <- mean((predictRF(er_servo,servo,checkCases=TRUE)-servo$Response)^2)
  mseAR_servo[sim] <- mean((predictRF(ar_servo,servo,checkCases=TRUE)-servo$Response)^2)
}

outData <- data.frame(dataset=rep("Servo", 4*nsim), sim=rep(1:nsim, 4), method=rep(c("RF", "SSS", "ER", "AR"), each=nsim),
                      mse = c(mseRF_servo, mseSSS_servo, mseER_servo, mseAR_servo))

#write.table(outData, file = "Results/predServo.csv", sep = ",", row.names=FALSE)





