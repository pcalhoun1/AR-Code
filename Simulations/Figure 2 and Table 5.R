#rm(list=ls(all=TRUE))


### Figure 2 ###
simulsNull <- read.csv("VIMP/Results/simulsNull.csv")

#setEPS()
#postscript(file="fig2.eps", width=5.07874, height=5.07874)
par(mfrow=c(2,2))

algoBoxplot <- function(data){
  par(mar=c(2,4,0,0)+0.1)
  rfBox<-boxplot(data, plot=FALSE)
  atX<-1:5
  bxp(rfBox,at=atX,boxwex=0.7,xlim=c(0.7,5.3),ylim=c(-0.08,0.08),axes=FALSE)
  axis(1,at=atX,labels=paste0("X",1:5),tick=FALSE,font=2,line=-1.5)
  axis(2,at=seq(-0.08,0.08,0.04),las=1,lwd=2,font=2)
  mtext("VIMP",2,font=2,line=2.8,cex=1.3)
}

algoBoxplot(simulsNull[,grep("rf",names(simulsNull))])
title("RF",line=-1.2,cex.main=2)

algoBoxplot(simulsNull[,grep("sss",names(simulsNull))])
title("SSS",line=-1.2,cex.main=2)

algoBoxplot(simulsNull[,grep("er",names(simulsNull))])
title("ER",line=-1.2,cex.main=2)

algoBoxplot(simulsNull[,grep("ar",names(simulsNull))])
title("AR",line=-1.2,cex.main=2)
#dev.off()


### Table 5 ###
propCorrect<-function(dataset){
  dat<-read.csv(paste0("VIMP/Results/",dataset))
  propCorrect_RF<-mean(apply(dat[,grep("rf",names(simulsNull))],1,function(x){which.max(x)==3}))
  propCorrect_SSS<-mean(apply(dat[,grep("sss",names(simulsNull))],1,function(x){which.max(x)==3}))
  propCorrect_ER<-mean(apply(dat[,grep("er",names(simulsNull))],1,function(x){which.max(x)==3}))
  propCorrect_AR<-mean(apply(dat[,grep("ar",names(simulsNull))],1,function(x){which.max(x)==3}))

  data.frame(propCorrect_RF,propCorrect_SSS,propCorrect_ER,propCorrect_AR)
}

data.frame(rel0.05=t(propCorrect("simuls0_05.csv")),
           rel0.10=t(propCorrect("simuls0_10.csv")),
           rel0.15=t(propCorrect("simuls0_15.csv")),
           rel0.20=t(propCorrect("simuls0_20.csv")))


