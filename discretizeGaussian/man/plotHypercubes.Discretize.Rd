%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotHypercubes.Discretize}
\alias{plotHypercubes.Discretize}
\alias{Discretize.plotHypercubes}
\alias{plotHypercubes,Discretize-method}

\title{Plotting the discretization of a bivariate random normal variable  }

\description{
  Plotting the discretization of a bivariate random normal variable 
}

\usage{
\method{plotHypercubes}{Discretize}(this, cubes, text="center", borders=FALSE, colors=NULL, ...)
}

\arguments{
 \item{cubes}{ The list of hypercubes. }
 \item{text}{ Text to be added to each hypercube. Value 'center'
   show the center point, 'index' show the index of the cube and if text i a vector of same length as the number of cube plot the text. }
 \item{borders}{ Show the border of the hypercubes if true.}
 \item{colors}{ A integer vector of same length as the number of cubes used to give the cubes colors. The color is set by the integer value. }
 \item{...}{Not used.}
}

\value{
  A plot is produced.
}

\author{Lars Relund \email{lars@relund.dk}}

\seealso{
  For more information see \code{\link{Discretize}}.
}

\examples{
dis<-Discretize()

## univariate
mu<-0
sigma2<-3
n<-10
dis$discretize1DUnifEqLth(mu,sigma2,n)  # eq length intervals
dis$discretize1DUnifEqProb(mu,sigma2,n) # equal probability

## bivariate (nonuniform)
mu<-c(0,0)
sigma<-matrix(c(2,-2,-2,4),nrow=2)

par(mfrow=c(2,3))
cubes<-dis$discretize2DNonunif(mu,sigma)
dis$plotHypercubes(cubes)
dis$plotHypercubes(cubes,borders=TRUE)
cols<-1:length(cubes)
dis$plotHypercubes(cubes, colors=cols)
cubes<-dis$discretize2DNonunif(mu,sigma,maxKL=0.3)
dis$plotHypercubes(cubes)
dis$plotHypercubes(cubes, text="index")
text<-paste("i",1:length(cubes)-1)
dis$plotHypercubes(cubes, text=text)

par(mfrow=c(2,2))
sigma<-matrix(c(4,0,0,4),nrow=2)
cubes<-dis$discretize2DNonunif(mu,sigma)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DNonunif(mu,sigma, modifyCenter="split1")
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DNonunif(mu,sigma, modifyCenter="split2")
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DNonunif(mu,sigma, modifyCenter="split2", split=0.5)
dis$plotHypercubes(cubes)

# split cubes further
par(mfrow=c(2,2))
cubes<-dis$discretize2DNonunif(mu,sigma)
dis$plotHypercubes(cubes, text="index")
cubes<-dis$splitCube2D(cubes,mu,sigma,iM=4)
dis$plotHypercubes(cubes, text="index")
cubes<-dis$splitCube2D(cubes,mu,sigma,iM=19)
dis$plotHypercubes(cubes, text="index")
cubes<-dis$splitCube2D(cubes,mu,sigma,iM=13)
dis$plotHypercubes(cubes, text="index")


## bivariate (uniform)
mu<-c(0,0)
sigma<-matrix(c(4,0,0,4),nrow=2)

par(mfrow=c(2,2))
cubes<-dis$discretize2DUnifEqInv(mu,sigma,4,6)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DUnifEqInv(mu,sigma,10,10)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DUnifEqProb(mu,sigma,4,6)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DUnifEqProb(mu,sigma,10,10)
dis$plotHypercubes(cubes)






}
\keyword{internal}
\keyword{methods}
