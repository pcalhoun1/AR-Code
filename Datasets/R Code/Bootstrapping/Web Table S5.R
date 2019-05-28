#rm(list=ls(all=TRUE))

source('../../../R Functions/RF functions 20JAN19.R')
load(file="../../Data/contData.RData")

### Change the RF functions to also return the number of attempts to split the node using the AR algorithm ###

pickCutpt <- function(allVars, y, ids, method, minsplit, minbucket, nGuesses, minpvalue, corstr){
  
  # Note: it is possible no candidate splits satisfy p-value threshold.  Thus, code will try nGuesses (default 50) before stopping
  for (guess in 1:nGuesses) {
    v <- sample.int(NCOL(allVars), 1)
    #If variable selected is a factor, pick a random subset of categories.
    #Do not sort by mean.  Otherwise, will be very likely to select categorical variable
    xTemp <- ordinalize(x=allVars[,v], y, ids, sortCat=FALSE, corstr=corstr)
    x <- xTemp$x
    #If all x values the same, do not check optimal split
    if (abs(max(x) - min(x)) > 1e-8) {
      cutpts <- findCutpts(x, minbucket)
      if (!is.null(cutpts)) {
        randCutpt <- ifelse(length(cutpts) == 1, cutpts, sample(cutpts,1))
        #Check whether random cutpoint satisfies quality control
        cutpt <- pvalue_check(y=y, x=x, ids=ids, cutpt=randCutpt, method=method, minpvalue=minpvalue, corstr=corstr)
        if (!is.null(cutpt)) {
          if (is.factor(allVars[,v])) {return(list(varid=v, cutpt=cutpt, x=x, cutToLvl=xTemp$cutToLvl, nTries=guess))
          } else {return(list(varid=v, cutpt=cutpt, x=x, cutToLvl=xTemp$cutToLvl, nTries=guess))}
        }
      }
    }
  }
}


partitionAR <- function(allVars, y, ids, subset, search, method, split, minsplit, minbucket, minpvalue, corstr){
  if (search != "ar") {stop("search not AR?")}
  if (sum(subset) < minsplit) {return(NULL)}
  allVars <- allVars[subset,,drop=FALSE]
  y <- y[subset]
  ids <- ids[subset]
  
  if (NROW(allVars) < 2*minbucket) {stop("Can't split tree to satisfy minbucket")}
  #If all y values are the same, rpart will give an error message.  Do not split tree
  if (length(unique(y)) == 1) {return(NULL)}
  
  randPt <- pickCutpt(allVars, y, ids, method, minsplit, minbucket, nGuesses=50, minpvalue=minpvalue, corstr=corstr)
  
  #It is possible (unlikely) no cutpoint selected or none can satisfy threshold criterion
  stat <- NA; cutoff <- NA
  if (!is.null(randPt)) {
    mod <- splitrule(y=y, x=randPt$x, ids=ids, cutpts=randPt$cutpt, method=method, split=split, corstr=corstr)
    stat <- mod$stat
    if(is.factor(allVars[,randPt$varid])){
      breakLeft <- rep(NA, length(levels(allVars[,randPt$varid])))
      breakLeft[levels(allVars[,randPt$varid]) %in% colnames(randPt$cutToLvl)[randPt$cutToLvl <= mod$cutoff]]=1L
      breakLeft[levels(allVars[,randPt$varid]) %in% colnames(randPt$cutToLvl)[randPt$cutToLvl > mod$cutoff]]=2L
      if(all(is.na(breakLeft)) & length(unique(breakLeft))<=1){stop("Did not find correct cutpoints")}
    } else {cutoff <- mod$cutoff}
  }
  
  #If each candidate variable cannot be split (e.g. cannot satisfy minbucket), return null
  if (is.na(stat)) {return(sum(subset))}
  if (is.na(cutoff)) {
    #Index is used for categorical variable splits
    return(partysplit(varid=randPt$varid,
                      index=breakLeft,
                      info=list(stats=stat, nTries=randPt$nTries)))
  } else {
    #Breaks is used for continuous variable splits
    return(partysplit(varid=randPt$varid,
                      breaks=cutoff,
                      info=list(stats=stat, nTries=randPt$nTries)))
  }
}


### ALL: Grow tree by using partition() function several times in a recursive loop ###
growTemp <- function(id=1L, depth=1L, data, response, idVar, subset, search, method, split,
                     mtry, nsplit, nsplit.random, minsplit, minbucket, maxdepth,
                     a, scale.y, useSearch, useOptim, useRpart, minpvalue, corstr){
  #print(c(depth,id))
  if (depth > maxdepth) {return(partynode(id=id))}
  
  y <- data[[response]]
  ids <- data[[idVar]]
  
  if (search=="ar") {
    sp <- partitionAR(allVars=data[,!(names(data) %in% c(response,idVar)),drop=FALSE], y=y, ids=ids, subset=subset, search=search,
                      method=method, split=split, minsplit=minsplit, minbucket=minbucket, minpvalue=minpvalue, corstr=corstr)
  } else {
    #Select candidate variables
    varSelected <- sort(sample.int(ncol(data)-1-(idVar != ""), mtry))
    vars <- data[varSelected]
    colnames(vars) <- varSelected #Have columns represent varid
    
    sp <- partition(vars=vars, y=y, ids=ids, subset=subset,
                    search=search, method=method, split=split, nsplit=nsplit, nsplit.random=nsplit.random,
                    minsplit=minsplit, minbucket=minbucket, a=a, scale.y=scale.y,
                    useSearch=useSearch, useOptim=useOptim, useRpart=useRpart, corstr=corstr)
  }
  
  if (is.null(sp)) {return(partynode(id=id))}
  if (length(sp)==1) {return(partynode(id=id, info=list(nTries=50, nObsAt50=sp)))}
  
  # Split the data
  kidids <- kidids_split(sp, data=data)
  depth <- depth + 1
  
  kids <- vector(mode="list", length=max(kidids, na.rm=TRUE))
  for (kidid in seq_along(kids)) {
    s <- subset
    s[kidids != kidid] <- FALSE
    # Node ID
    if (kidid > 1) {myid <- max(nodeids(kids[[kidid-1]]))
    } else {myid <- id}
    # Start recursion on this daugther node
    kids[[kidid]] <- growTemp(id=as.integer(myid+1), depth=depth, data=data, response=response, idVar=idVar,
                              subset=s, search=search, method=method, split=split, mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random,
                              minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                              a=a, scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart,
                              minpvalue=minpvalue, corstr=corstr)
  }
  return(partynode(id=as.integer(id), split=sp, kids=kids,
                   info=list(stats=min(info_split(sp)$stats, na.rm=TRUE), nTries=info_split(sp)$nTries)))
}

# Function to count the depth of a node of a tree #
library(stringr)
idDepth <- function(tree) {
  outTree <- capture.output(tree)
  idCount <- 1
  depthValues <- rep(NA, length(tree))
  names(depthValues) <- 1:length(tree)
  for (index in seq_along(outTree)){
    if (grepl("\\[[0-9]+\\]", outTree[index])) {
      depthValues[idCount] <- str_count(outTree[index], "\\|")
      idCount = idCount + 1
    }
  }
  return(depthValues)
}




### Friedman1 ###

ntrees <- 100

# friedman1 #
friedman1 <- contData$friedman1
form <- as.formula(paste0("Response ~ ",paste0("v",1:10,collapse=' + ')))

set.seed(486194)
ar_friedman1 <- growRF(ntrees=ntrees, formula=form, data=friedman1, search="ar", method="anova", split="MSE",
                       mtry=1, minsplit=6, minbucket=3, maxdepth=10, minpvalue=0.05, sampleMethod='bootstrap')

# Track number of attempts made before splitting the node
# Also track terminal nodes where the minsplit was satisfied, but couldn't find split after 50 guesses
nTriesInternal <- NULL
nTriesTerminal <- NULL
for (treeNum in 1:ntrees) {
  ar_friedman1Tree <- ar_friedman1[[treeNum]]$tree
  
  # Count number of guesses were tried for each split
  nTriesInt <- data.frame(internalNodes = nodeids(ar_friedman1Tree)[!(nodeids(ar_friedman1Tree) %in% nodeids(ar_friedman1Tree, terminal = TRUE))])
  depthAll <- idDepth(ar_friedman1[[treeNum]]$tree)
  nTriesInt$depth <- depthAll[names(depthAll) %in% nTriesInt$internalNodes]
  nTriesInt$nTries <- unlist(nodeapply(ar_friedman1Tree, ids = nTriesInt$internalNodes, FUN = function(n) n$info$nTries))
  nTriesInternal <- rbind(nTriesInternal, nTriesInt)
  
  # Count number of times 50 guesses were tried and no split was chosen.
  # For these cases, look at number of observations and depth of terminal node
  nObsTerm <- unlist(nodeapply(ar_friedman1Tree, ids = nodeids(ar_friedman1Tree)[(nodeids(ar_friedman1Tree) %in% nodeids(ar_friedman1Tree, terminal = TRUE))], FUN = function(n) n$info$nObsAt50))
  depthTerm <- depthAll[names(nObsTerm)]
  nTriesTerminal <- rbind(nTriesTerminal, data.frame(id=names(nObsTerm), nObsTerm=nObsTerm, depthTerm=depthTerm))
  rownames(nTriesTerminal) <- NULL
}



#head(nTriesDat)
nTriesGrp <- rep("", nrow(nTriesInternal))
nTriesGrp[1 <= nTriesInternal$nTries & nTriesInternal$nTries <= 2] <- "A: 1-2"
nTriesGrp[3 <= nTriesInternal$nTries & nTriesInternal$nTries <= 5] <- "B: 3-5"
nTriesGrp[6 <= nTriesInternal$nTries & nTriesInternal$nTries <= 10] <- "C: 6-10"
nTriesGrp[11 <= nTriesInternal$nTries & nTriesInternal$nTries <= 20] <- "D: 11-20"
nTriesGrp[21 <= nTriesInternal$nTries & nTriesInternal$nTries <= 30] <- "E: 21-30"
nTriesGrp[31 <= nTriesInternal$nTries & nTriesInternal$nTries <= 40] <- "F: 31-40"
nTriesGrp[41 <= nTriesInternal$nTries & nTriesInternal$nTries <= 50] <- "G: 41-50"
nTriesGrp <- c(nTriesGrp, rep("H: No Split", nrow(nTriesTerminal)))
depthValues <- c(nTriesInternal$depth, nTriesTerminal$depth) + 1

table(depthValues)

### Table S5 ###
round(prop.table(table(depthValues, nTriesGrp), 1)*100, 1)

median(nTriesTerminal$nObsTerm)  # Median number of observations with no split
mean(nTriesTerminal$nObsTerm)  # Mean number of observations with no split
max(nTriesTerminal$nObsTerm)  # Max number of observations with no split


