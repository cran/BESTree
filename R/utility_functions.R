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
## This file contains all function that should not be called by users directly.
## These functions are used by callable functions.
###################################################################################

# ImpFct provides the vector Impurity (compile only one of the following)
# Input : a vector of categorical variable
# Output : Impurity value for the vector
ImpFctCM <- function(CM) {



  #return(((CM[,2])/length(Y))%*%(1-(CM[,2])/length(Y)))

  return(1-sum(((CM[,2])/sum(CM[,2]))^2))


}

#ImpFct <- compiler::cmpfun(ImpFctCM)

#Deviance (without count)
#ImpFctD <- function(CM) {

#  return( -(log((CM[,2])/length(Y)) %*% (CM[,2])/length(Y)) )


#}

# BestNum finds the Best Split within numerical predictors
# BestNum2 Is an attempt at removing counting from loops
# Input : Data (Data frame, last column is response)
# Output : Information about best split among numerical predictors (Vector)
BestNum <- function(Data) {


  n <- nrow(Data)
  d <- ncol(Data)-1

  X <- data.frame(Data[,1:d])
  Y <- as.integer(factor(Data[,(d+1)]))

  CM <- plyr::count(Y)

  Impurity <- n*ImpFctCM(CM)
  NImpurity <- 0

  #Value to be returned
  SV <- 0
  SP <- 0
  ImpRed <- 0

  for ( j in 1:(d)) {

    #Ordering is really heavy computation wise :(
    Data <- Data[order(Data[,j]),]
    X <- data.frame(Data[,1:d])
    Y <- as.integer(factor(Data[,(d+1)]))

    XV <- Data[1,j]
    Index <- 1

    CM1 <- matrix(data=rep(0,nrow(CM)*2),ncol=2)
    CM1[as.integer(Y[1]),2] <- 1

    CM2 <- CM

    CM2[as.integer(Y[1]),2] <- CM2[as.integer(Y[1]),2] - 1

    while ( Index < n-6) {

      if ( round(Data[(Index+1),j],10) != round(XV,10) && Index > 5 ) {

        NImpurity <- (Index*ImpFctCM(CM1)+(n-Index)*ImpFctCM(CM2))



        if ( NImpurity < Impurity ) {

          SV <- j
          SP <- (XV +  Data[(Index+1),j])/2
          ImpRed <- ImpRed + (Impurity - NImpurity)
          ind <- Index
          Impurity <- NImpurity

        }

        XV <- Data[(Index+1),j]


      }

      Index <- Index+1

      CM1[as.integer(Y[Index]),2] <- CM1[as.integer(Y[Index]),2] +1

      CM2[as.integer(Y[Index]),2] <- CM2[as.integer(Y[Index]),2] - 1
    }
  }


  BNum <- c(SV,SP,ImpRed)

  return(BNum)


}



# OrderOnce will contain a matrix, containing the ordered index for numerical predictors
# Input : Data is a data frame containing the numerical predictors
# Ouput : Matrix? containg order based upon every numerical predictors
OrderOnce <- function(Data) {

  d <- ncol(Data)
  n <- nrow(Data)

  Order <- matrix(data=rep(0,d*n),nrow=d,ncol=n)

  for ( j in 1:d ) {
    Order[j,] <- order(Data[,j])

  }

  return(Order)
}

# BestCat finds the Best Split within categorical predictors
# Input : Data (Data frame, last column is response)
# Output : Information about best split among categorical predictors (List)
BestCat <- function(Data,RegionIndex,CatList, CatVANo) {

  n <- nrow(Data[RegionIndex,])
  d <- ncol(Data)-1

  X <- data.frame(Data[,1:d])
  Y <- Data[,(d+1)]

  CM <- plyr::count(Data[RegionIndex,(d+1)])

  Impurity <- n*ImpFctCM(CM)
  ImpRed <- 0
  SV <- 0
  LCat <- 0

  for ( j in 1:d) {

    NoLvl <- length(levels(factor(X[RegionIndex,j])))
    NoSplit <- 2^(NoLvl-1)-1


    if ( NoSplit > 0 ) {

      #First split ( These set operators are SUPER SLOW)
      IndexL <- intersect( RegionIndex, CatList[[CatVANo[j]]][[as.integer(levels(factor(X[RegionIndex,j]))[NoLvl])]] )
      IndexR <- setdiff( RegionIndex, IndexL )

      NImpurity <- length(IndexL)*ImpFctCM(plyr::count(Data[IndexL,(d+1)]))+length(IndexR)*ImpFctCM(plyr::count(Data[IndexR,(d+1)]))

      if ( NImpurity < Impurity ) {

        SV <- j
        LCat <- levels(factor(X[RegionIndex,j]))[NoLvl]
        ImpRed <- ImpRed + (Impurity - NImpurity)

        Impurity <- NImpurity

      }

      if ( NoLvl > 2 ) {

        for ( s in 1:(NoLvl-2)) {

          ToAdd <- utils::combn((NoLvl-1),s)

          for ( i in 1:ncol(ToAdd)) {

            NIndexL <- IndexL

            # Add new category in Left leaf ( These set operators are SUPER SLOW)
            for ( k in 1:s) {
              NIndexL <- c(NIndexL,intersect( RegionIndex, CatList[[CatVANo[j]]][[ToAdd[k,i]]] ))
            }

            NIndexR <- setdiff( RegionIndex, NIndexL )

            NImpurity <- length(NIndexL)*ImpFctCM(plyr::count(Y[NIndexL]))+length(NIndexR)*ImpFctCM(plyr::count(Y[NIndexR]))

            if ( NImpurity < Impurity ) {

              SV <- j
              LCat <- c(NoLvl,ToAdd[,i])
              ImpRed <- ImpRed + (Impurity - NImpurity)

              Impurity <- NImpurity

            }

          }


        }
      }

    }


  }

  ToReturn <- list()
  ToReturn[[1]] <- SV
  ToReturn[[2]] <- LCat
  ToReturn[[3]] <- ImpRed

  return(ToReturn)
}

# CatIndex returns a list containing index for categorical values
# Input : Data is a data frame containing the categorical predictors
# Output : List where first increment is categorical variable number
# and 2nd sub list is level value
CatIndex <- function(Data,CatNo) {

  d <- ncol(Data)
  n <- nrow(Data)

  Index <- list()

  for ( j in 1:d ) {
    Index[[CatNo[j]]] <- list()
    for ( i in 1:length(levels(Data[,j])) ) {
      Index[[CatNo[j]]][[i]] <- (Data[,j]==levels(Data[,j])[i])*(1:n)
      Index[[CatNo[j]]][[i]] <-  Index[[CatNo[j]]][[i]][ !is.na(Index[[CatNo[j]]][[i]]) ]
      Index[[CatNo[j]]][[i]] <-  Index[[CatNo[j]]][[i]][ Index[[CatNo[j]]][[i]] >0 ]
    }

  }

  return(Index)
}

# ListOfTree Builds all possible pruned trees that will be used in the tree pruning process
# Input : A BEST Object
# Output : a list of Trees (Matrix form)
ListOfTrees <- function(Fit){


  Tree <- data.frame(Fit[[1]])
  Regions <- Fit[[3]]

  TreeList <- list()

  TreeList[[1]] <- Tree

  NoLeaves <- sum(Tree[,9]==1)

  while( NoLeaves > 1 ) {

    PossibleC <- Tree[Tree[,9]==0,]

    CInd <- rep(0,nrow(PossibleC))

    for ( i in 1:nrow(PossibleC)) {

      if ( Tree[PossibleC[i,4],9] ==1 && Tree[PossibleC[i,5],9] ==1 ) {

        CInd[i] <- 1
      }


    }

    RPC <- PossibleC[CInd==1,]

    ToC <- RPC[which.min(RPC[,8]),]

    #2 stands for 'collapsed leaves'
    Tree[ToC[1,4],9] <- 2
    Tree[ToC[1,5],9] <- 2

    Tree[ToC[1,1],9] <- 1

    TreeList[[length(TreeList)+1]] <- Tree

    NoLeaves <- NoLeaves-1


  }

  return(TreeList)
}
