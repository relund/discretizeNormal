%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Discretize}
\docType{class}
\alias{Discretize}


\title{Discretize class}

\description{
 Containing all methods related to discritizing an uni/bi-variate normal distribution. Both uniform and nonuniform discretization possible.
 Package:  discretizeGaussian \cr
\bold{Class Discretize}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Discretize}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Discretize}\cr
extends \link[R.oo]{Object}\cr


}

\usage{
Discretize(...)
}

\arguments{
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[discretizeGaussian:addCube.Discretize]{addCube}} \tab Adds a 2D cube to the plot.\cr
 \tab \code{\link[discretizeGaussian:addCubeCol.Discretize]{addCubeCol}} \tab Adds a 2D cube with color to the plot.\cr
 \tab \code{\link[discretizeGaussian:addIdx.Discretize]{addIdx}} \tab Add cube index to the plot.\cr
 \tab \code{\link[discretizeGaussian:addPoints.Discretize]{addPoints}} \tab Adds center points to the plot.\cr
 \tab \code{\link[discretizeGaussian:addText.Discretize]{addText}} \tab Add text to the plot.\cr
 \tab \code{\link[discretizeGaussian:bounds1D.Discretize]{bounds1D}} \tab Min, mean and max density on a 1D cube.\cr
 \tab \code{\link[discretizeGaussian:bounds2D.Discretize]{bounds2D}} \tab Min, mean and max density on a 2D cube.\cr
 \tab \code{\link[discretizeGaussian:direc.Discretize]{direc}} \tab Finds the optimal (approximate) direcection to spilt a cube.\cr
 \tab \code{\link[discretizeGaussian:discretize1DUnifEqLth.Discretize]{discretize1DUnifEqLth}} \tab Discretize a normal distribution such that intervals have equal length.\cr
 \tab \code{\link[discretizeGaussian:discretize1DUnifEqProb.Discretize]{discretize1DUnifEqProb}} \tab Discretize a normal distribution such that intervals have equal probability.\cr
 \tab \code{\link[discretizeGaussian:discretize1DVec.Discretize]{discretize1DVec}} \tab Discretize the real numbers according to a set of center points.\cr
 \tab \code{\link[discretizeGaussian:discretize2DNonunif.Discretize]{discretize2DNonunif}} \tab Discretize a bivariate normal distribution using a non-uniform discretization .\cr
 \tab \code{\link[discretizeGaussian:discretize2DUnifEqInv.Discretize]{discretize2DUnifEqInv}} \tab Discretize a bivariate normal distribution using a uniform discretization with intervals of equal length .\cr
 \tab \code{\link[discretizeGaussian:discretize2DUnifEqProb.Discretize]{discretize2DUnifEqProb}} \tab Discretize a bivariate normal distribution using a uniform discretization with intervals of equal probability .\cr
 \tab \code{\link[discretizeGaussian:klBound1D.Discretize]{klBound1D}} \tab Upper bound on KL distance on a 1D cube.\cr
 \tab \code{\link[discretizeGaussian:klBound2D.Discretize]{klBound2D}} \tab Upper bound on KL distance on a 2D cube.\cr
 \tab \code{\link[discretizeGaussian:plotCubes.Discretize]{plotCubes}} \tab Plot the cubes (only bivariate distributions).\cr
 \tab \code{\link[discretizeGaussian:plotHypercubes.Discretize]{plotHypercubes}} \tab Plotting the discretization of a bivariate random normal variable  .\cr
 \tab \code{\link[discretizeGaussian:ratio.Discretize]{ratio}} \tab Calc max divided by min density value.\cr
 \tab \code{\link[discretizeGaussian:splitCube2D.Discretize]{splitCube2D}} \tab Split a cube further up .\cr
 \tab \code{\link[discretizeGaussian:volCube.Discretize]{volCube}} \tab Volume/length of cube.\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


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

\references{
[1] Nielsen, L.R.; Jørgensen, E. & Højsgaard, S. Embedding a state space model into a Markov decision process Dept. of Genetics and Biotechnology, Aarhus University, 2008. \cr
}

\author{Lars Relund \email{lars@relund.dk}}
\keyword{classes}
