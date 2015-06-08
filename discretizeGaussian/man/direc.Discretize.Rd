%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{direc.Discretize}
\alias{direc.Discretize}
\alias{Discretize.direc}
\alias{direc,Discretize-method}

\title{Finds the optimal (approximate) direcection to spilt a cube}

\description{
  Internal function.
}

\usage{
\method{direc}{Discretize}(this, cube, mu, sigma, ...)
}

\arguments{
 \item{cube}{ The cube under consideration which is a (2x2) matrix containing the bounds of the variables (columnwise) (bivariate case).   }
 \item{mu}{ The mean. }
 \item{sigma}{ The covariate matrix. }
 \item{...}{Not used.}
}

\value{
  Return the variable index to split.
}

\author{Lars Relund \email{lars@relund.dk}}

\seealso{
  For more information see \code{\link{Discretize}}.
}


\keyword{internal}
\keyword{methods}
