#rm(list=ls(all=TRUE))

### PARALLEL PROGRAMMING ###

library(parallel)

getwd()
dir()

singleSim<-function(x, ntrees, relevance){

  source('../../R Functions/RF functions 20JAN19.R')  

  x1 <- rnorm(120,mean=0,sd=1)
  x2 <- as.factor(sample(letters[1:2], 120, replace=TRUE, prob=rep(0.5,2)))
  x3 <- as.factor(sample(letters[1:4], 120, replace=TRUE, prob=rep(0.5,4)))
  x4 <- as.factor(sample(letters[1:10], 120, replace=TRUE, prob=rep(0.5,10)))
  x5 <- as.factor(sample(letters[1:20], 120, replace=TRUE, prob=rep(0.5,20)))
  
  Class <- rep(NA,120)
  Class[x3 %in% c("a","b")] <- rbinom(sum(x3 %in% c("a","b")),size=1,prob=.5-relevance)
  Class[x3 %in% c("c","d")] <- rbinom(sum(x3 %in% c("c","d")),size=1,prob=.5+relevance)
  Class <- as.factor(Class)
  
  dat <- data.frame(x1,x2,x3,x4,x5,Class)
  form <- as.formula("Class ~ x1 + x2 + x3 + x4 + x5")
  
  rf_vimp <- vimpRF(growRF(ntrees=ntrees, formula=form, data=dat, search="exhaustive", method="class", split="gini",
                            mtry=2, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE),
                    dat,checkCases=TRUE)$scaledVI
  sss_vimp <- vimpRF(growRF(ntrees=ntrees, formula=form, data=dat, search="sss", method="class", split="gini",
                             mtry=2, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap'),
                     dat,checkCases=TRUE)$scaledVI
  er_vimp <- vimpRF(growRF(ntrees=ntrees, formula=form, data=dat, search="exhaustive", method="class", split="gini",
                            mtry=2, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap'),
                    dat,checkCases=TRUE)$scaledVI
  ar_vimp <- vimpRF(growRF(ntrees=ntrees, formula=form, data=dat, search="ar", method="class", split="gini",
                            mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=0.05, sampleMethod='bootstrap'),
                    dat,checkCases=TRUE)$scaledVI
  
  return(list(rf_vimp=rf_vimp,sss_vimp=sss_vimp,er_vimp=er_vimp,ar_vimp=ar_vimp))
}


simVimpPar<-function(nsim, ntrees, relevance){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim")
  clusterSetRNGStream(cl=cl, iseed=684602)
  simuls <- parLapply(cl,1:nsim,singleSim,ntrees=ntrees,relevance=relevance)
  stopCluster(cl)

  rf_vimp <- matrix(NA,ncol=5,nrow=nsim)
  sss_vimp <- matrix(NA,ncol=5,nrow=nsim)
  er_vimp <- matrix(NA,ncol=5,nrow=nsim)
  ar_vimp <- matrix(NA,ncol=5,nrow=nsim)
  for(sim in 1:nsim){
    rf_vimp[sim,] <- simuls[[sim]]$rf_vimp
    sss_vimp[sim,] <- simuls[[sim]]$sss_vimp
    er_vimp[sim,] <- simuls[[sim]]$er_vimp
    ar_vimp[sim,] <- simuls[[sim]]$ar_vimp
  }
  return(list(rf_vimp=rf_vimp, sss_vimp=sss_vimp, er_vimp=er_vimp, ar_vimp=ar_vimp))
}

simulsNull <- simVimpPar(nsim=1000, ntrees=2000, relevance=0)
#write.table(simulsNull, file = "Results/simulsNull.csv", sep = ",")

simuls0_05 <- simVimpPar(nsim=1000, ntrees=2000, relevance=0.05)
#write.table(simuls0_05, file = "Results/simuls0_05.csv", sep = ",")

simuls0_10 <- simVimpPar(nsim=1000, ntrees=2000, relevance=0.10)
#write.table(simuls0_10, file = "Results/simuls0_10.csv", sep = ",")

simuls0_15 <- simVimpPar(nsim=1000, ntrees=2000, relevance=0.15)
#write.table(simuls0_15, file = "Results/simuls0_15.csv", sep = ",")

simuls0_20 <- simVimpPar(nsim=1000, ntrees=2000, relevance=0.20)
#write.table(simuls0_20, file = "Results/simuls0_20.csv", sep = ",")

