#rm(list=ls(all=TRUE))

library(dplyr)
library(reshape2)

##### BOOTSTRAP #####

fileList <- list.files("Vimp and Prediction/Results")

results <- NULL
for (file in fileList){
  results <-  bind_rows(results, read.csv(paste0("Vimp and Prediction/Results/", file)))
}
dim(results)

head(results)

results$model[results$model == "Interaction"] <- "1.Interaction"
results$model[results$model == "Binary Split"] <- "2.Binary Split"
results$model[results$model == "Nonlinear"] <- "3.Nonlinear"

results$type[results$type == "quantitative"] <- "1.Quantitative"
results$type[results$type == "categorical"] <- "2.Categorical"

# 10 categorical predictors shows basically the same trends as 10 quantitative predictors.  Don't report for simplicity
results <- results[!(results$type=="2.Categorical" & results$nValues==10), ]

# Prediction #
aucResults <- aggregate(list(results[["AUC"]]), list(results[["model"]], results[["nvars"]], results[["nValues"]], results[["type"]],  results[["input"]], results[["method"]]),
                         FUN=function(x){c(mean=round(mean(x),3), sd=round(sd(x),6))})
aucResults <- data.frame(model=aucResults[,1], nvars=aucResults[,2], nValues=aucResults[,3], type=aucResults[,4], method=aucResults[,5], input=aucResults[,6],
                         aucMean=aucResults[,7][,1], aucSD = aucResults[,7][,2])

# VIMP #

results$rankX1 <- apply(results, 1, function(x) {rank(-as.numeric(x[grepl("x", names(x))]), ties.method = "first")[1]})
results$rankX2 <- apply(results, 1, function(x) {rank(-as.numeric(x[grepl("x", names(x))]), ties.method = "first")[2]})
results$top2 <- apply(results, 1, function(x) {max(as.numeric(c(x[["rankX1"]], x[["rankX2"]]))) <= 2})

vimpResults <- aggregate(list(results[["top2"]]), list(results[["model"]], results[["nvars"]], results[["type"]], results[["nValues"]], results[["input"]], results[["method"]]),
                         FUN=function(x){c(pctTop2 = round(mean(x)*100,0), nTop2 = sum(x))})
vimpResults <- data.frame(model=vimpResults[,1], nvars=vimpResults[,2], type=vimpResults[,3], nValues=vimpResults[,4], method=vimpResults[,5], input=vimpResults[,6],
                          pctTop2=vimpResults[,7][,1], nTop2 = vimpResults[,7][,2])

### Table 7 ###
dcast(aucResults, model + nvars + type + nValues ~ method + input, value.var="aucMean")[c(1:4, 7,8,6,5)]  # Default - AUC
dcast(vimpResults, model + nvars + type + nValues ~ method + input, value.var="pctTop2")[c(1:4, 7,8,6,5)]  # Default - VIMP


