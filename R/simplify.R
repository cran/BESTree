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
## This file contains two speciall callable functions.
## These are meant to make your life easier.
###################################################################################

#' Computes the proportion of matching terms in two vectors of the same length.
#' Used to compute the accuracy for prediction on test set.
#' @param Vec1 A vector of labels
#' @param Vec2 Another vector of labels
#' @return Percentage of identical labels (accuracy)
#' @examples
#' Vec1 <- c(1,1,2,3,1)
#' Vec2 <- c(1,2,2,3,1)
#' Acc(Vec1,Vec2)
#' @export
Acc <- function(Vec1,Vec2) {

  Transformed <- c(Vec1,Vec2)

  Vec1I <- as.integer(Transformed[1:length(Vec1)])
  Vec2I <- as.integer(Transformed[(length(Vec1)+1):(2*length(Vec1))])

  acc <- sum((Vec1I-Vec2I)==0)/length(Vec1)

  return(acc)

}

#' Quickly build the Available Variable list necessary for BEST
#' This list contains details as to which variables is available for the partitioning.
#' It also contains which variables are gating variables.
#' @param d  Number of predictors
#' @param GV Gating variables
#' @param BEV Branch-Exclusive Variables
#' @param Thresh Threshold for Gates
#' @param Direc Direction of Gates ( 1 means add variable if bigger than thresh)
#' @return The list containing the Variable Availability structure
#' @examples
#' #This function can be used to set up the variable availability structure.
#' #Suppose we want to fit a regular decision tree on a data set containing d predictors
#' d <- 10
#' VA <- ForgeVA(d,1,0,0,0)
#' #Suppose now that predictor x5 is a binary gating variable for x4
#' #such that x4 is available if x5 = 1
#' GV <- 5 #The gating variable
#' BEV <- 4 #The Branch-Exclusive variable
#' Tresh = 0.5 #Value between 0 and 1
#' Direc = 1 #X4 is available if X5 is bigger than Tresh
#' VA <- ForgeVA(d,GV,BEV,Tresh,Direc)
#' @export
ForgeVA <- function(d,GV,BEV,Thresh=0.5,Direc=0) {

  VA <- list()

  VA[[1]] <- rep(1,d)

  VA[[1]][BEV] <- 0

  for ( i in 1:d) {

    VA[[i+1]] <- list()

    VA[[i+1]][[1]] <- 0

    VA[[i+1]][[2]] <- rep(0,d)

    VA[[i+1]][[3]] <- rep(0,d)


  }

  for ( j in 1:length(GV) ) {

    VA[[GV[j]+1]][[1]] <- Thresh[j]

    if ( Direc[j] == 1 ) {

      VA[[GV[j]+1]][[3]][BEV[j]] <- 1
    } else {
      VA[[GV[j]+1]][[2]][BEV[j]] <- 1
    }



  }

  return(VA)
}
