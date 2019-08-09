###################################################################################
###   Branch-Exclusive Splits Trees
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
## This file contains callable functions related to simplest version of the BEST
## algorithm.
###################################################################################

###################################################################################
## BEST
##
## Parameters
## Data - Data set (Data Frame) : Can take on both numerical and categorical
## predictors. Last Column is Response Variable (Categorical only)
## Integer needed for factor levels
## Size - Minimal Number of Observation within a leaf needed for partitionning
## VA - Variable Availability structure
##
## Returns
## Tree - Matrix representing a Tree including Region Names, No Obs, Split
## Variables, Split points, prediction for each region, etc...
## Regions - Row Names for every Regions
## Split Points - Necessar for categorical predictors
###################################################################################

#' Main function of the package.
#' It produces Classification Trees with Branch-Exclusive variables.
#' @param Data A data set (Data Frame): Can take on both numerical and categorical predictors. Last column of the data set must be the Repsonse Variable (Categorical Variables only)
#' @param Size Minimal Number of Observation within a leaf needed for partitionning
#' @param VA Variable Availability structure
#' @return A BEST object with is a list containing the resulting tree, row numbers for each regions and the split points
#' @examples
#' n <- 1000
#' Data <- BESTree::Data[1:n,]
#' d <- ncol(Data)-1
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' Fit <- BESTree::BEST(Data,Size,VA)
#' @export
BEST <- function( Data, Size, VA ) {

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

    Response <- plyr::count(Y)

    Prediction <- Response$x[which.max(Response$freq)]
    Impurity <- n*ImpFctCM(Response)

    # Preparing Variable Available
    CatVA <- Variables[[Ite]]*CatInd
    NumVA <- (Variables[[Ite]]-CatVA)

    CatVANo <- which(CatVA > 0)
    NumVANo <- which(NumVA > 0)

    # Making sure splitting IS possible
    if ( sum(CatVA)>0 && n > 1) {
      CatCond <- FALSE
      for ( i in 1:length(CatVANo)) {
        CatCond  <- CatCond  || length(unique(X[,CatVANo[i]])) > 1
      }


    } else {
      CatCond<- FALSE
    }

    if ( sum(NumVA)>0 && n > 1) {
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
        BN <- BestNum(data.frame(Data[Regions[[Ite]],c(NumVANo,(d+1))]))
      } else {
        BN <- c(0,0,0)
      }

      #Proceed to split on Categorical predictors if possible
      if ( CatCond ) {
        BC <- BestCat(data.frame(Data[,c(CatVANo,(d+1))]),Regions[[Ite]],CatList, CatVANo)
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

        Regions[[LB]] <- intersect( Regions[[Ite]], unlist(CatList[[SV]][as.integer(SP)]) )

        Regions[[BB]] <- setdiff(Regions[[Ite]],Regions[[LB]])

        Tree <- rbind(Tree,c(Ite,SV,-1,LB,BB,Prediction,n,BC[[3]],0,Impurity) )

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

#' Classify a new observation point
#' @param Point A new observation
#' @param Fit A BEST object
#' @return The predicted class
#' @examples
#' n <- 500
#' Data <- BESTree::Data[1:n,]
#' NewPoint <- BESTree::Data[n+1,]
#' d <- ncol(Data)-1
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' Fit <- BESTree::BEST(Data,Size,VA)
#' BESTree::Predict(NewPoint[1:d],Fit)
#' @export
Predict <- function(Point,Fit) {

  Tree <- Fit[[1]]
  SPList <- Fit[[2]]
  CatInd <- Fit[[4]]
  Reg <- 1
  Prediction <- 0
  LN <- Tree[Reg,9]

  while ( LN == 0 ) {

    SV <- Tree[Reg,2]
    #Check if SV is Categorical
    if ( CatInd[SV] == 1 ) {

      if ( is.element(Point[SV],SPList[[Reg]]) ) {

        Reg <- Tree[Reg,4]

      } else {

        Reg <- Tree[Reg,5]

      }


    # Else the SV is numerical
    } else {

      if ( Point[SV] <= Tree[Reg,3] ) {
        Reg <- Tree[Reg,4]
      } else {
        Reg <- Tree[Reg,5]
      }
    }

    LN <- Tree[Reg,9]

  }


  Prediction <- Tree[Reg,6]

  return(Prediction)
}

#' Classify a set of new observation points
#' @param M A matrix of new observations where one row is one observation
#' @param Fit A BEST object
#' @return The predicted class
#' @examples
#' n <- 500
#' Data <- BESTree::Data[1:n,]
#' d <- ncol(Data)-1
#' NewPoints <- BESTree::Data[(n+1):(n+11),1:d]
#' VA <- ForgeVA(d,1,0,0,0)
#' Size <- 50
#' Fit <- BESTree::BEST(Data,Size,VA)
#' Predictions <- BESTree::MPredict(NewPoints,Fit)
#' @export
MPredict <- function(M,Fit) {

  n <- nrow(M)

  Predictions <- rep(0,n)

  for ( i in 1:n) {

    Predictions[i] <- Predict(M[i,],Fit)
  }



  return(Predictions)
}

#' Uses a Validation Set to select the best trees within the list of pruned trees.
#' @param Fit A BEST object
#' @param VSet A Validation Set (Can also be used in CV loop)
#' @return The shallower trees among trees wiht Highest accuracy. This replaces the first element in the BEST object list.
#' @examples
#' nv <- 50
#' ValData <- BESTree::Data[(1000+1):nv,]
#' Fit <- BESTree::Fit
#' Fit[[1]] <- BESTree::TreePruning(Fit,ValData)
#' @export
TreePruning <- function(Fit,VSet){

  #List of possible pruned trees
  TreeList <- ListOfTrees(Fit)

  #Number of trees to test
  NTree <- length(TreeList)

  d <- ncol(VSet)-1

  Accuracy <- rep(0,NTree)

  Pred <- MPredict(VSet[,(1:d)],Fit)

  Accuracy[1] <- Acc(VSet[,(d+1)],Pred)

  #Set accuracy of full tree as Benchmark
  BestAccuracy <- Accuracy[1]

  Index <- 1

  for ( i in 2:NTree) {

    Fit[[1]] <- TreeList[[i]]

    Pred <- MPredict(VSet[,(1:d)],Fit)

    Accuracy[i] <- Acc(VSet[,(d+1)],Pred)

    #If one tree beat benchmark then it is the new best tree
    #We use greater or EQUAL since we want trees as shallow as possible
    if ( Accuracy[i] >= BestAccuracy ) {

      BestAccuracy <- Accuracy[i]

      Index <- i

    }

  }



  return(TreeList[[Index]])


}

