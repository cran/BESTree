% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bestforest.R
\name{VI}
\alias{VI}
\title{Produces a variable important analysis using the mean decrease in node impurity}
\usage{
VI(Forest)
}
\arguments{
\item{Forest}{A list of BEST Objects (Usually produced by RBEST or BESTForest)}
}
\value{
A vector of importance (size d)
}
\description{
Produces a variable important analysis using the mean decrease in node impurity
}
\examples{
n <- 500
Data <- BESTree::Data[1:n,]
d <- ncol(Data)-1
NewPoints <- BESTree::Data[(n+1):(n+11),1:d]
VA <- ForgeVA(d,1,0,0,0)
Size <- 50
NoT <- 10
Fit <- BESTree::BaggedBEST(Data,VA,NoT,Size)
VI <- BESTree::VI(Fit)
}
