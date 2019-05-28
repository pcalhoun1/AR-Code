#rm(list=ls(all=TRUE))
library(mlbench)


#######################################################################
###                        BINARY DATASETS                          ###
###                                                                 ###
#######################################################################

### Breast cancer dataset ###
data(BreastCancer)
BreastCancer <- BreastCancer[complete.cases(BreastCancer),] #Remove 16 observations with a missing value

#Convert ordinal factors to numeric
BreastCancer$Cl.thickness <- as.numeric(BreastCancer$Cl.thickness)
BreastCancer$Cell.size <- as.numeric(BreastCancer$Cell.size)
BreastCancer$Cell.shape <- as.numeric(BreastCancer$Cell.shape)
BreastCancer$Marg.adhesion <- as.numeric(BreastCancer$Marg.adhesion)
BreastCancer$Epith.c.size <- as.numeric(BreastCancer$Epith.c.size)
BreastCancer$Class <- as.factor(BreastCancer$Class)
dim(BreastCancer); str(BreastCancer)

### Sonar dataset ###
data(Sonar)
dim(Sonar); str(Sonar)
any(is.na(Sonar))  #Note, no missing values

### Ionosphere dataset ###
data(Ionosphere)
dim(Ionosphere); str(Ionosphere)
any(is.na(Ionosphere))  #Note, no missing values

### German credit dataset ###
german <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data", 
                     header=FALSE, nrows=1000)
colnames(german)[21] <- "Class"
german$Class <- as.factor(german$Class)
dim(german); str(german)
any(is.na(german))  #Note, no missing values

### Liver dataset ###
liver <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/liver-disorders/bupa.data", 
                    header=FALSE, sep=",", nrows=345)
colnames(liver) <- c("mcv", "alkphos", "sgpt", "sgot", "gammagt", "drinks", "Class")
liver$Class <- as.factor(liver$Class)
dim(liver); str(liver)
any(is.na(liver))  #Note, no missing values

### Twonorm dataset ###
set.seed(871320)
temp <- mlbench.twonorm(3300)
twonorm <- data.frame(cbind(temp$x,temp$classes))
colnames(twonorm) <- c(paste0("V",1:20),"Class")
twonorm$Class <- as.factor(twonorm$Class)
dim(twonorm); str(twonorm)
any(is.na(twonorm))  #Note, no missing values

### Threenorm dataset ###
set.seed(413732)
temp <- mlbench.threenorm(3300)
threenorm <- data.frame(cbind(temp$x,temp$classes))
colnames(threenorm) <- c(paste0("V",1:20),"Class")
threenorm$Class <- as.factor(threenorm$Class)
dim(threenorm); str(threenorm)
any(is.na(threenorm))  #Note, no missing values

### Ringnorm dataset ###
set.seed(904536)
temp <- mlbench.ringnorm(3300)
ringnorm <- data.frame(cbind(temp$x,temp$classes))
colnames(ringnorm) <- c(paste0("V",1:20),"Class")
ringnorm$Class <- as.factor(ringnorm$Class)
dim(ringnorm); str(ringnorm)
any(is.na(ringnorm))  #Note, no missing values

### Spambase dataset ###
spambase <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.data", 
                       header=FALSE, sep=",")
colnames(spambase)[58] <- "Class"
spambase$Class <- as.factor(spambase$Class)
dim(spambase); str(spambase)
any(is.na(spambase))  #Note, no missing values

### Birthwt dataset ###
data(birthwt, package="MASS")
birthwt <- birthwt[,1:9] #Remove bwt variable
colnames(birthwt)[1] <- "Class"
birthwt$Class <- as.factor(birthwt$Class)
birthwt$race <- as.factor(birthwt$race)
birthwt$smoke <- as.factor(birthwt$smoke)
birthwt$ht <- as.factor(birthwt$ht)
birthwt$ui <- as.factor(birthwt$ui)
dim(birthwt); str(birthwt)
any(is.na(birthwt))  #Note, no missing values


### Combine all binary datasets ###
binData <- list(BreastCancer, Sonar, Ionosphere, german, liver, twonorm, threenorm, ringnorm, spambase, birthwt)
names(binData) <- c("breastCancer", "sonar", "ionosphere", "german", "liver", "twonorm", "threenorm", "ringnorm", "spambase", "birthwt")

#save(binData, file="binData.RData")
#load(file="binData.RData")


#######################################################################
###                      Regression DATASETS                        ###
###                                                                 ###
#######################################################################

### Housing dataset ###
data(BostonHousing2)
BostonHousing2 <- BostonHousing2[c(1:4,7:19,6)]  #cmedv is response (remove medv), rearrange dataset
names(BostonHousing2)[18] <- "Response"
dim(BostonHousing2); str(BostonHousing2)
any(is.na(BostonHousing2)) #Note, no missing values

### Servo dataset ###
data(Servo)
names(Servo)[5] <- "Response"
dim(Servo); str(Servo)
any(is.na(Servo))  #Note, no missing values but only factors

### Abalone dataset ###
abalone <- read.table('http://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data', header = FALSE , sep = ',',
                     col.names = c('sex','length','diameter','height','weight.w','weight.s','weight.v','weight.sh','response'))
names(abalone)[9] <- "Response"
dim(abalone); str(abalone)
any(is.na(abalone)) #Note, no missing values

### Friedman1 dataset ###
set.seed(308501)
temp <- mlbench.friedman1(2200)
friedman1 <- data.frame(cbind(temp$x,temp$y))
colnames(friedman1) <- c(paste0("v",1:10),"Response")
dim(friedman1); str(friedman1)
any(is.na(friedman1)) #Note, no missing values

### Friedman2 dataset ###
set.seed(768108)
temp <- mlbench.friedman2(2200)
friedman2 <- data.frame(cbind(temp$x,temp$y))
colnames(friedman2) <- c(paste0("v",1:4),"Response")
dim(friedman2); str(friedman2)
any(is.na(friedman2)) #Note, no missing values

### Friedman3 dataset ###
set.seed(279661)
temp <- mlbench.friedman3(2200)
friedman3 <- data.frame(cbind(temp$x,temp$y))
colnames(friedman3) <- c(paste0("v",1:4),"Response")
dim(friedman3); str(friedman3)
any(is.na(friedman3)) #Note, no missing values

### Ailerons dataset ###
ailerons <- read.csv("ailerons.csv",header=FALSE)
colnames(ailerons)[41]="Response"
dim(ailerons); str(ailerons)
any(is.na(ailerons)) #No missing values

### Elevators dataset ###
elevatorsTrain <- read.table("elevators.data", header = FALSE , sep = ',')
elevatorsTest <- read.table("elevators.test", header = FALSE , sep = ',')
elevators <- rbind(elevatorsTrain,elevatorsTest)
colnames(elevators)[19] <- "Response"
dim(elevators); str(elevators)
any(is.na(elevators)) #No missing values

### imports85 ###
data(imports85, package="randomForest")
imports85 <- imports85[,-2] # Too many NAs in normalizedLosses.
imports85 <- imports85[complete.cases(imports85), ]  #Remove missing data
colnames(imports85)[25] <- "Response"
#Convert ordinal factors to numeric
imports85$numOfCylinders <- as.numeric(imports85$numOfCylinders)
dim(imports85); str(imports85)
any(is.na(imports85)) #No missing values

### data(airquality) ###
data(airquality)
airquality <- airquality[complete.cases(airquality), ]  #Remove missing data
colnames(airquality)[1]="Response"
dim(airquality); str(airquality)
any(is.na(airquality)) #No missing values


### Combine all continuous datasets ###
contData <- list(BostonHousing2, Servo, abalone, friedman1, friedman2, friedman3, ailerons, elevators, imports85, airquality)
names(contData) <- c("housing", "servo", "abalone", "friedman1", "friedman2", "friedman3", 'ailerons', "elevators", "imports85", "airquality")

#save(contData, file="contData.RData")
#load(file="contData.RData")



