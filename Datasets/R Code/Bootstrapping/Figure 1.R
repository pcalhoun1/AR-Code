#rm(list=ls(all=TRUE))
getwd()
dir()

source('../../../R Functions/RF functions 20JAN19.R')
load(file="../../Data/binData.RData")

# Sonar dataset with RF and 1,000 trees takes around 2 minutes; manuscript repeats this 50 times #
ntrees=1000

# Sonar #
sonar <- binData$sonar
form <- as.formula(paste0("Class ~ ",paste(paste0("V",1:60),collapse=' + ')))
# Note: Table 2 calculates AUC 50 times and takes the average.  To demonstrate what the AUC looks like,
# pick a seed with an AUC very close to the average across 50 times (AUC = 0.9345) #
rf_sonar <- growRF_Parallel(ntrees=ntrees, formula=form, data=sonar, search="exhaustive", method="class", split="gini",
                            mtry=8, nsplit=NULL, minsplit=6, minbucket=3, maxdepth=30, sampleMethod='bootstrap', useRpart=TRUE,
                            iseed=121)
lvls <- levels(sonar$Class)
predRF <- predictRF(rf_sonar,sonar,checkCases=TRUE)
# Show AUC for this seed #
auc(sonar$Class, predRF, direction="<")[1]

#setEPS()
#postscript(file="fig1.eps", width=2.75591, height=2.75591, pointsize=8)
par(mar=c(3, 3.2, 0, 0.1) + 0.1)
plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
mtext("1-Specificity", 1, line=1.8, cex=1.7, font=2)
mtext("Sensitivity", 2, line=2, cex=1.7, font=2)
axis(1, seq(0, 1, .2), line=-0.4, font=2, mgp=c(3, .7, 0), cex.axis=1.3)
axis(2, seq(0, 1, .2),las=2, line=-0.3, font=2, mgp=c(3, .7, 0), cex.axis=1.3)
xCurve <- roc(sonar$Class, predRF, direction="<")$specificities
yCurve <- roc(sonar$Class, predRF, direction="<")$sensitivities
lines(x=1-xCurve, y=yCurve, lwd=2)
polygon(c(0,1,1-xCurve),c(0,0,yCurve),density=7, col='black')
legend(0.5, 0.35, "AUC = 0.9345", pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')
dev.off()
