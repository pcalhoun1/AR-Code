#rm(list=ls(all=TRUE))

### PARALLEL PROGRAMMING ###

library(parallel)

getwd()
dir()

singleSim<-function(x, model = "Interaction", ntrees=100, nUselessPred = 3, type = "quantitative", nValues = 10){
  
  source('../../R Functions/RF functions 20JAN19.R')  
  
  x1 <- sample(c(-5:-1, 1:5), 1100, replace=TRUE)
  x2 <- sample(c(-5:-1, 1:5), 1100, replace=TRUE)

  if (model == "Interaction") {linPred <- 0 + 2*x1 + 2*x2 - 4*x1*x2
  } else if (model == "Binary Split") {linPred <- 0 + 3*(x1 <= -3) + 3*(x2 <= -3) - 6*(x1 <= -3)*(x2 <= -3) - 3*(x1 >= 3) - 3*(x2 >= 3) + 6*(x1 >= 3)*(x2 >= 3)
  } else if (model == "Nonlinear") {linPred <- 2*sin(pi/8*x1) + 2*sin(pi/8*x2)
  } else {stop("Unexpected Model")}
  
  prob = 1/(1+exp(-linPred))
  Class <- rbinom(1100, 1, prob)
  Class <- as.factor(Class)
  #table(Class)
  if (type == "quantitative") { dat <- data.frame(Class, x1, x2)
  } else if (type == "categorical") { dat <- data.frame(Class, x1=as.factor(x1), x2=as.factor(x2))
  } else { stop("Unexpected type?") }
  
  for (p in 3:(2+nUselessPred)) {
    dat[[paste0("x",p)]] <- sample((-nValues/2):(nValues/2-1), 1100, replace=TRUE)
    if (type == "categorical") { dat[[paste0("x",p)]] <- as.factor(dat[[paste0("x",p)]]) }
  }
  
  form <- as.formula(paste0("Class ~ ", paste0("x", 1:(2+nUselessPred), collapse=" + ")))

  aucRF <- aucSSS <- aucER <- rep(NA, (2+nUselessPred))
  for (mtry in 1:(2+nUselessPred)) {
    # Split data into training and test dataset #
    assign(paste0("rf_sim", mtry),
           growRF(ntrees=ntrees, formula=form, data=dat[1:100, ], search="exhaustive", method="class", split="gini",
                     mtry=mtry, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE))
    assign(paste0("sss_sim", mtry),
           growRF(ntrees=ntrees, formula=form, data=dat[1:100, ], search="sss", method="class", split="gini",
                      mtry=mtry, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, a=50, sampleMethod='bootstrap'))
    assign(paste0("er_sim", mtry),
           growRF(ntrees=ntrees, formula=form, data=dat[1:100, ], search="exhaustive", method="class", split="gini",
                     mtry=mtry, nsplit=1, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap'))
    lvls <- levels(dat$Classdat[1:100])
    
    predRF <- predictRF(get(paste0("rf_sim", mtry)), dat[1:100, ], checkCases=TRUE)
    predSSS <- predictRF(get(paste0("sss_sim", mtry)), dat[1:100, ], checkCases=TRUE)
    predER <- predictRF(get(paste0("er_sim", mtry)), dat[1:100, ], checkCases=TRUE)
    
    aucRF[mtry] <- auc(dat$Class[1:100], predRF, direction="<")[1]
    aucSSS[mtry] <- auc(dat$Class[1:100], predSSS, direction="<")[1]
    aucER[mtry] <- auc(dat$Class[1:100], predER, direction="<")[1]
  }
  
  outAUC <- NULL
  outVimp <- NULL
  
  mtry_default <- round(sqrt(2+nUselessPred))
  for (method in c("RF", "SSS", "ER")) {
    mtry_optim <- which.max(get(paste0("auc", method)))
    pred_default <- predictRF(get(paste0(tolower(method), "_sim", mtry_default)), dat[101:1100, ], checkCases=TRUE)
    pred_optim <- predictRF(get(paste0(tolower(method), "_sim", mtry_optim)), dat[101:1100, ], checkCases=TRUE)
    auc_default <- auc(dat$Class[101:1100], pred_default, direction="<")[1]
    auc_optim <- auc(dat$Class[101:1100], pred_optim, direction="<")[1]
    outAUC <- c(outAUC, auc_default, auc_optim)
    
    # Use training dataset to get vimp #
    vimp_default <- vimpRF(get(paste0(tolower(method), "_sim", mtry_default)), dat[1:100, ], checkCases=TRUE)$scaledVI
    vimp_optim <- vimpRF(get(paste0(tolower(method), "_sim", mtry_optim)), dat[1:100, ], checkCases=TRUE)$scaledVI
    outVimp <- rbind(outVimp, vimp_default, vimp_optim)
  }

  
  aucAR <- rep(NA, 6)
  index <- 1
  minpvalues <- c(0.01, seq(0.05, 0.25, 0.05))
  for (minpvalue in minpvalues) {
    # Split data into training and test dataset #
    assign(paste0("ar_sim", minpvalue),
           growRF(ntrees=ntrees, formula=form, data=dat[1:100, ], search="ar", method="class", split="gini",
                  mtry=1, minsplit=6, minbucket=3, maxdepth=30, minpvalue=minpvalue, sampleMethod='bootstrap'))

    lvls <- levels(dat$Classdat[1:100])
    predAR <- predictRF(get(paste0("ar_sim", minpvalue)), dat[1:100, ], checkCases=TRUE)
    aucAR[index] <- auc(dat$Class[1:100], predAR, direction="<")[1]
    index <- index + 1
  }
  
  minpvalue_optim <- which.max(aucAR)
  pred_default <- predictRF(ar_sim0.05, dat[101:1100, ], checkCases=TRUE)
  pred_optim <- predictRF(get(paste0("ar_sim", minpvalues[minpvalue_optim])), dat[101:1100, ], checkCases=TRUE)
  auc_default <- auc(dat$Class[101:1100], pred_default, direction="<")[1]
  auc_optim <- auc(dat$Class[101:1100], pred_optim, direction="<")[1]
  outAUC <- c(outAUC, auc_default, auc_optim)
  
  # Use training dataset to get vimp #
  vimp_default <- vimpRF(ar_sim0.05, dat[1:100, ], checkCases=TRUE)$scaledVI
  vimp_optim <- vimpRF(get(paste0("ar_sim", minpvalues[minpvalue_optim])), dat[1:100, ], checkCases=TRUE)$scaledVI
  outVimp <- rbind(outVimp, vimp_default, vimp_optim)
  rownames(outVimp) <- NULL

  return(data.frame(model=model, ntrees = ntrees, nvars = 2 + nUselessPred, type=type, nValues = nValues,
                    method=rep(c("RF","SSS","ER","AR"),each=2), input = rep(c("default", "optim"), 4), AUC=outAUC, outVimp))
}


simVimpPar<-function(nsim, model, ntrees, nUselessPred, type, nValues){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim")
  clusterSetRNGStream(cl=cl, iseed=684602)
  simuls <- parLapply(cl,1:nsim, singleSim, model=model, ntrees=ntrees, nUselessPred=nUselessPred, type=type, nValues=nValues)
  stopCluster(cl)
  
  out <- data.frame(sim = rep(1:nsim, each=8))
  out <- cbind(out, do.call("rbind",simuls))
  
  return(out)
}

simulsIntA <- simVimpPar(nsim=500, model="Interaction", ntrees=100, nUselessPred=3, type="quantitative", nValues=10)
#write.table(simulsIntA, file = "Results/simulsIntA.csv", sep = ",", row.names=FALSE)

simulsIntB <- simVimpPar(nsim=500, model="Interaction", ntrees=100, nUselessPred=8, type="quantitative", nValues=10)
#write.table(simulsIntB, file = "Results/simulsIntB.csv", sep = ",", row.names=FALSE)

simulsIntC <- simVimpPar(nsim=500, model="Interaction", ntrees=100, nUselessPred=8, type="quantitative", nValues=50)
#write.table(simulsIntC, file = "Results/simulsIntC.csv", sep = ",", row.names=FALSE)

simulsIntD <- simVimpPar(nsim=500, model="Interaction", ntrees=100, nUselessPred=8, type="categorical", nValues=50)
#write.table(simulsIntD, file = "Results/simulsIntD.csv", sep = ",", row.names=FALSE)

simulsIntE <- simVimpPar(nsim=500, model="Interaction", ntrees=100, nUselessPred=8, type="categorical", nValues=10)
#write.table(simulsIntE, file = "Results/simulsIntE.csv", sep = ",", row.names=FALSE)



simulsBinA <- simVimpPar(nsim=500, model="Binary Split", ntrees=100, nUselessPred=3, type="quantitative", nValues=10)
#write.table(simulsBinA, file = "Results/simulsBinA.csv", sep = ",", row.names=FALSE)

simulsBinB <- simVimpPar(nsim=500, model="Binary Split", ntrees=100, nUselessPred=8, type="quantitative", nValues=10)
#write.table(simulsBinB, file = "Results/simulsBinB.csv", sep = ",", row.names=FALSE)

simulsBinC <- simVimpPar(nsim=500, model="Binary Split", ntrees=100, nUselessPred=8, type="quantitative", nValues=50)
#write.table(simulsBinC, file = "Results/simulsBinC.csv", sep = ",", row.names=FALSE)

simulsBinD <- simVimpPar(nsim=500, model="Binary Split", ntrees=100, nUselessPred=8, type="categorical", nValues=50)
#write.table(simulsBinD, file = "Results/simulsBinD.csv", sep = ",", row.names=FALSE)

simulsBinE <- simVimpPar(nsim=500, model="Binary Split", ntrees=100, nUselessPred=8, type="categorical", nValues=10)
#write.table(simulsBinE, file = "Results/simulsBinE.csv", sep = ",", row.names=FALSE)



simulsNonlinearA <- simVimpPar(nsim=500, model="Nonlinear", ntrees=100, nUselessPred=3, type="quantitative", nValues=10)
#write.table(simulsNonlinearA, file = "Results/simulsNonlinearA.csv", sep = ",", row.names=FALSE)

simulsNonlinearB <- simVimpPar(nsim=500, model="Nonlinear", ntrees=100, nUselessPred=8, type="quantitative", nValues=10)
#write.table(simulsNonlinearB, file = "Results/simulsNonlinearB.csv", sep = ",", row.names=FALSE)

simulsNonlinearC <- simVimpPar(nsim=500, model="Nonlinear", ntrees=100, nUselessPred=8, type="quantitative", nValues=50)
#write.table(simulsNonlinearC, file = "Results/simulsNonlinearC.csv", sep = ",", row.names=FALSE)

simulsNonlinearD <- simVimpPar(nsim=500, model="Nonlinear", ntrees=100, nUselessPred=8, type="categorical", nValues=50)
#write.table(simulsNonlinearD, file = "Results/simulsNonlinearD.csv", sep = ",", row.names=FALSE)

simulsNonlinearE <- simVimpPar(nsim=500, model="Nonlinear", ntrees=100, nUselessPred=8, type="categorical", nValues=10)
#write.table(simulsNonlinearE, file = "Results/simulsNonlinearE.csv", sep = ",", row.names=FALSE)
