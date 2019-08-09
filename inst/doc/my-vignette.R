## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(BESTree)

## ----include = TRUE------------------------------------------------------
set.seed(100)
n=1000
X1 <- rnorm(n,0,sd=1)
X2 <- rnorm(n,2,sd=2)
X3 <- runif(n,0,1)
X4 <- runif(n,-2,2)

Y <- 1*(X1<0)*(X4<0.5)+0*(X1>0)*(X4<0.5)+1*(X3>0.5)*(X4>0.5)+0*(X3<0.5)*(X4>0.5)
#Add some randomized Y
RY <- sample(1000,150)
Y[RY] <- 1-Y[RY]

## ----include = TRUE------------------------------------------------------

X3[X3>0.5] <- NA
Data <- cbind(X1,X2,X3,X4,as.factor(Y))


## ----include = TRUE------------------------------------------------------

X5 <- is.na(X3)*1
NewData <- cbind(Data[,1:4],X5,Data[,ncol(Data)])

Training <- NewData[1:800,]
Valid <- NewData[801:900,]
Testing <- NewData[901:1000,]

d = ncol(NewData)-1 #number of predictor
VA <- BESTree::ForgeVA(d,5,3)

## ----include = TRUE------------------------------------------------------
VA

## ----include = TRUE------------------------------------------------------
Fit <- BESTree::BEST(Training,10,VA)
PTree <- BESTree::TreePruning(Fit,Valid)
Fit[[1]] <- PTree
preds <- BESTree::MPredict(Testing[,1:d],Fit)
BESTree::Acc(preds,Testing[,d+1])

