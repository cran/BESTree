###################################################################################
###   Branch-Exclusive Splits Trees (Non-Callables)
###################################################################################
###   Cedric Beaulac (2018)
###   Fully Re-Coded BEST Algorithm (Classification Tree only)
###################################################################################
###   Can only manage simple gating structure ( January 23rd, 2018 )
###################################################################################
###   Bagged Tree, Random Forest and Variable Importances  ( February 27th, 2018 )
###################################################################################
###   R Packaging  ( April 5th, 2019 )
###################################################################################

###################################################################################
## This file contains all callable functions related to Forest and Tree Bagging
###################################################################################

#' Performs Bootstrap Aggregating of BEST trees
#' @param Data A data set (Data Frame): Can take on both numerical and categorical predictors. Last column of the data set must be the Repsonse Variable (Categorical Variables only)
#' @param VA Variable Availability structure
#' @param Size Minimal Number of Observation within a leaf needed for partitionning (default is 50)
#' @param NoT Number of Trees in the bag
#' @return A list of BEST Objects
#' @examples
#' n <- 500
#' Data <- BESTree::Data[1:n,]
#' d <- ncol(Data)-1
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' NoT <- 10
#' Fit <- BESTree::BaggedBEST(Data,VA,NoT,Size)
#' @export
BaggedBEST <- function(Data,VA,NoT=50,Size=50) {

  #To be returned is a list of BEST objects
  ListofFit <- list()
  n <- nrow(Data)

  for ( i in 1:NoT ) {

    #Bootstrap sample
    samp <- sample((1:n),size=n,replace=TRUE)

    BData <- Data[samp,]

    #Fit BEST Trees, 50 is abitrary
    BFit <- BEST(BData,Size,VA)

    ListofFit[[i]] <- BFit


  }


  return(ListofFit)
}

# Construct a BEST tree where predictors are randomly selected (for RF purposes)
# Input : A training Set, the size of trees and Variable availability structure
# Ouptut : A BEST Objects
RBEST <- function( Data, Size, VA ) {

  # Establishing basic data information
  n <- nrow(Data)
  d <- (ncol(Data) - 1)

  X <- data.frame(Data[,1:d])
  Y <- Data[,(d+1)]

  #Determine if a input is factor or not
  CatInd <- rep(0,d)
  for ( i in 1:d ) {
    if (is.factor(X[,i])) {
      CatInd[i] <- 1
    }
  }

  #Set up list of index for categorical predictors
  CatNo <- CatInd*(1:d)
  if ( sum(CatNo) > 0 ) {
    CatList <- CatIndex(data.frame(Data[,CatNo]),setdiff(CatNo,0))
  }


  # Innteger indication the current iteration
  Ite <- 1

  # Cond is a binary variable to check if we continue splitting or not
  Cond <- TRUE

  #Variable avaible at any iteration
  Variables <- list()
  Variables[[Ite]] <- VA[[1]]

  #Regions
  Regions <- list()
  Regions[[Ite]] <- seq(1:n)

  #Tree object
  Tree <- matrix( data = c(0,0,0,0,0,0,0,0,1,1), nrow=1 )

  #Split Points, this is mostly for Categorical variables purposes, since we can't include
  #vectors in the Tree matrix
  SPList <- list()

  while ( Ite <= length(Regions) ) {

    n <- length(Regions[[Ite]])


    X <- data.frame(Data[Regions[[Ite]],1:d])
    Y <- Data[Regions[[Ite]],(d+1)]

    Response <-  plyr::count(Y)

    Prediction <- Response$x[which.max(Response$freq)]
    Impurity <- n*ImpFctCM(Response)

    #Randomly picking variables ( For random forest purposes)
    dIte <- sum(Variables[[Ite]])
    ToDraw <- round(sqrt(dIte))
    From <- which(Variables[[Ite]] > 0)

    Var <- sample(From,size=ToDraw,replace=FALSE)



    # Preparing Variable Available
    CatVA <- Variables[[Ite]]*CatInd
    NumVA <- (Variables[[Ite]]-CatVA)

    CatVANo <- which(CatVA > 0)
    NumVANo <- which(NumVA > 0)

    CatVANo <- intersect(CatVANo,Var)
    NumVANo <- intersect(NumVANo,Var)



    # Making sure splitting IS possible
    if ( length(CatVANo)>0 && n > 1) {
      CatCond <- FALSE
      for ( i in 1:length(CatVANo)) {
        CatCond  <- CatCond  || length(unique(X[,CatVANo[i]])) > 1
      }


    } else {
      CatCond<- FALSE
    }

    if ( length(NumVANo)>0 && n > 1) {
      NumCond <- FALSE
      for ( i in 1:length(NumVANo)) {
        NumCond  <- NumCond  || length(unique(X[,NumVANo[i]])) > 1
      }


    } else {
      NumCond <- FALSE
    }




    # Check Conditions
    Cond <- ( (n > Size) && nrow(Response) > 1 && (NumCond || CatCond))

    # If one of the conditions is false, we stop the splitting process on that branch.
    if ( !Cond ) {

      Tree <- rbind(Tree , c(Ite,0,0,0,0,Prediction,n,0,1,Impurity))

      SPList[[Ite]] <- 0

      #Tree[[Ite]] <- c(Ite,0,0,0,0,Pred,n,0,1)

    } else {

      #Proceed to split on Numerical predictors if possible
      if ( NumCond ) {
        BN <- BestNum(Data[Regions[[Ite]],c(NumVANo,(d+1))])
      } else {
        BN <- c(0,0,0)
      }

      #Proceed to split on Categorical predictors if possible
      if ( CatCond ) {
        BC <- BestCat(data.frame(Data[,c(CatVANo,(d+1))]),Regions[[Ite]],CatList,CatVANo)
      } else {
        BC <- list(0,0,0)
      }


      #Select the type of split producing best decreasin in impurity

      #If Numerical split produce best impurity reduction
      if ( BN[3] == 0 && BC[[3]] == 0 ) {

        Tree <- rbind(Tree, c(Ite,0,0,0,0,Prediction,n,0,1,Impurity))

        #Tree[[Ite]] <- c(Ite,0,0,0,0,Pred,n,0,1)

        SPList[[Ite]] <- 0

      } else if (BN[3] >= BC[[3]] ) {


        #Keep track of Split variable and point select
        #Also, names of children leaves and Impurity drop produced
        SV <- NumVANo[BN[1]]

        SP <- BN[2]

        LB <- length(Regions)+1

        BB <- length(Regions)+2

        Regions[[LB]] <- intersect(which(Data[,SV]<=SP),Regions[[Ite]])

        Regions[[BB]] <- setdiff(Regions[[Ite]],Regions[[LB]])

        Tree <- rbind(Tree,c(Ite,SV,SP,LB,BB,Prediction,n,BN[3],0,Impurity)  )

        #Tree[[Ite]] <- c(Ite,SV,SP,LB,BB,Pred,n,BN[3],0)

        SPList[[Ite]] <- SP

        #BEST : Update the list of predictors available for futur splits
        #If SP if greater then threshold only top data get new variables
        if ( SP > VA[[SV+1]][[1]] ) {


          Variables[[LB]] <- Variables[[Ite]]

          Variables[[BB]] <- Variables[[Ite]]+VA[[SV+1]][[3]]

          #If SP is smaller than threshold only lower data get new variables
        } else if (SP < VA[[SV+1]][[1]] ) {


          Variables[[LB]] <- Variables[[Ite]]+VA[[SV+1]][[2]]

          Variables[[BB]] <- Variables[[Ite]]

        } else {

          Variables[[LB]] <- Variables[[Ite]]+VA[[SV+1]][[2]]

          Variables[[BB]] <- Variables[[Ite]]+VA[[SV+1]][[3]]

        }


      } else if (BC[[3]] > BN[3] ) {

        #Keep track of Split variable and point select
        #Also, names of children leaves and Impurity drop produced
        SV <- CatVANo[BC[[1]]]

        SP <- BC[[2]]

        LB <- length(Regions)+1

        BB <- length(Regions)+2

        Regions[[LB]] <- intersect( Regions[[Ite]], unlist(CatList[[ SV ]][as.integer(SP)]) )

        Regions[[BB]] <- setdiff(Regions[[Ite]],Regions[[LB]])

        Tree <- rbind(Tree,c(Ite,SV,-1,LB,BB,Prediction,n,BN[[3]],0,Impurity) )

        #Tree[[Ite]] <- c(Ite,SV,-1,LB,BB,Pred,n,BN[[3]],0)

        #BEST : Update the list of predictors available for futur splits
        #Not coded yet, Gating variable MUST BE numerical

        Variables[[LB]] <- Variables[[Ite]]

        Variables[[BB]] <- Variables[[Ite]]

        SPList[[Ite]] <- SP

      }
    }


    Ite <- Ite+1

  }

  colnames(Tree) <- c('RegNo','SV','SP','LLeaf','RLeaf','Pred','NoObs','ImpRed','Leaf/Node','Imp')

  ToReturn <- list()

  ToReturn[[1]] <- Tree[2:nrow(Tree),]
  ToReturn[[2]] <- SPList
  ToReturn[[3]] <- Regions
  ToReturn[[4]] <- CatInd


  return(ToReturn)
}

#' Generates a random forest of BEST trees
#' @param Data A data set (Data Frame): Can take on both numerical and categorical predictors. Last column of the data set must be the Repsonse Variable (Categorical Variables only)
#' @param VA Variable Availability structure
#' @param Size Minimal Number of Observation within a leaf needed for partitionning (default is 50)
#' @param NoT Number of Trees in the bag
#' @return A list of BEST Objects (Random Forest)
#' @examples
#' n <- 500
#' Data <- BESTree::Data[1:n,]
#' d <- ncol(Data)-1
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' NoT <- 10
#' Fit <- BESTree::BESTForest(Data,VA,NoT,Size)
#' @export
BESTForest <- function (Data,VA,NoT=50,Size=50) {

  ListofFit <- list()
  n <- nrow(Data)

  for ( i in 1:NoT ) {

    samp <- sample((1:n),size=n,replace=TRUE)

    BData <- Data[samp,]

    BFit <- RBEST(BData,Size,VA)

    ListofFit[[i]] <- BFit


  }


  return(ListofFit)
}

#' Emits prediction from a forest of BEST's
#' @param M A matrix of new observations where one row is one observation
#' @param  LFit A list of BEST Objects (Usually produced by RBEST or BESTForest)
#' @return A vector of predictions
#' @examples
#' n <- 500
#' Data <- BESTree::Data[1:n,]
#' d <- ncol(Data)-1
#' NewPoints <- BESTree::Data[(n+1):(n+11),1:d]
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' NoT <- 10
#' Fit <- BESTree::BaggedBEST(Data,VA,NoT,Size)
#' Predictions <- BESTree::FPredict(NewPoints,Fit)
#' @export
FPredict <- function(M,LFit)  {

  NoT <- length(LFit)

  n <- nrow(M)

  Predictions <- rep(0,n)

  Pred <- matrix(rep(0,n*NoT),nrow=n)


  for ( i in 1:n) {

    for ( j in 1:NoT) {

      #Matrix of predictions base of every trees
      Pred[i,j] <- Predict(M[i,],LFit[[j]])
    }

    #Aggregating the predictions (Majority of votes)
    Predictions[i] <- plyr::count(Pred[i,])$x[which.max(plyr::count(Pred[i,])$freq)]

  }





  return(Predictions)
}

#' Produces a variable important analysis using the mean decrease in node impurity
#' @param Forest A list of BEST Objects (Usually produced by RBEST or BESTForest)
#' @return A vector of importance (size d)
#' @examples
#' n <- 500
#' Data <- BESTree::Data[1:n,]
#' d <- ncol(Data)-1
#' NewPoints <- BESTree::Data[(n+1):(n+11),1:d]
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' NoT <- 10
#' Fit <- BESTree::BaggedBEST(Data,VA,NoT,Size)
#' VI <- BESTree::VI(Fit)
#' @export
VI <- function(Forest){

  d <- length(Forest[[1]][[4]])

  NT <- length(Forest)

  TImp <- rep(0,d)


  for ( i in 1:NT) {

    ImpT <- rep(0,d)

    Tree <- Forest[[i]][[1]]

    DTree <- data.frame(Tree)

    #ddplyRes <- plyr::ddply(DTree,"SV",plyr::summarize,MeanImpRed = sum(DTree$ImpRed))

    ddplyRes <- stats::aggregate(DTree, list(SV = DTree$SV), sum)[,c(1,9)]

    for ( j in 2:nrow(ddplyRes)){

      ImpT[ddplyRes[j,1]] <- ddplyRes[j,2]
    }

    TImp <- TImp+ImpT
  }


  Imp <- TImp/NT

  return(Imp)
}
