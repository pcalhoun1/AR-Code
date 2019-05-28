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

mseAR_servo <- rep(NA, 6*nsim)
index <- 1
for (sim in 1:nsim) {
  for (minpvalue in c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90)) {
    ar_servo<-growRF_Parallel(ntrees=ntrees, formula=form, data=servo, search="ar", method="anova", split="MSE",
                              mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap', iseed=sim)
    mseAR_servo[index] <- mean((predictRF(ar_servo,servo,checkCases=TRUE)-servo$Response)^2)
    index <- index + 1
  }
}

outData <- data.frame(dataset = rep("Servo", 6*nsim), sim = rep(1:nsim, each=6), method = rep("AR", 6*nsim),
                      minpvalue = rep(c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90), nsim),
                      mse = mseAR_servo)

#write.table(outData, file = "Results/minpvaluePredServo.csv", sep = ",", row.names=FALSE)





