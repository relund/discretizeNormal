%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{discretize2DNonunif.Discretize}
\alias{discretize2DNonunif.Discretize}
\alias{Discretize.discretize2DNonunif}
\alias{discretize2DNonunif,Discretize-method}

\title{Discretize a bivariate normal distribution using a non-uniform discretization }

\description{
  Discretize a bivariate normal distribution into hypercubes (squares)
  such that the approximation have a certain Kulback Libler (KL) distance.
}

\usage{
\method{discretize2DNonunif}{Discretize}(this, mu, sigma, maxKL=0.5, maxIte=500, modifyCenter="no",
  split=0.25, ...)
}

\arguments{
 \item{mu}{ The mean (2-dim vector). }
 \item{sigma}{ The covariance (2x2 matrix). }
 \item{maxKL}{ Max KL distance. }
 \item{maxIte}{ Max number of iterations. }
 \item{modifyCenter}{ If no don't split the cubes around the mean center. If "split1" split the 4 cubes around the mean into 9 squares such that the mean is the center of a cube. If "split2" first add cubes such that the axis of the mean always in the center of the cubes. }
 \item{split}{ Only used if modifyCenter = "split2" to set the size of the nine cubes around the mean. }
 \item{...}{Not used.}
}

\value{
  A list where each element describe the cube and contains: KL - an upper bound on the KL-distance, cube - the bounds, center - the center, idxM - the index, cubeB - the fixed bounds (used for plotting).
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
