% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{ForgeVA}
\alias{ForgeVA}
\title{Quickly build the Available Variable list necessary for BEST
This list contains details as to which variables is available for the partitioning.
It also contains which variables are gating variables.}
\usage{
ForgeVA(d, GV, BEV, Thresh = 0.5, Direc = 0)
}
\arguments{
\item{d}{Number of predictors}

\item{GV}{Gating variables}

\item{BEV}{Branch-Exclusive Variables}

\item{Thresh}{Threshold for Gates}

\item{Direc}{Direction of Gates ( 1 means add variable if bigger than thresh)}
}
\value{
The list containing the Variable Availability structure
}
\description{
Quickly build the Available Variable list necessary for BEST
This list contains details as to which variables is available for the partitioning.
It also contains which variables are gating variables.
}
\examples{
#This function can be used to set up the variable availability structure.
#Suppose we want to fit a regular decision tree on a data set containing d predictors
d <- 10
VA <- ForgeVA(d,1,0,0,0)
#Suppose now that predictor x5 is a binary gating variable for x4
#such that x4 is available if x5 = 1
GV <- 5 #The gating variable
BEV <- 4 #The Branch-Exclusive variable
Tresh = 0.5 #Value between 0 and 1
Direc = 1 #X4 is available if X5 is bigger than Tresh
VA <- ForgeVA(d,GV,BEV,Tresh,Direc)
}
