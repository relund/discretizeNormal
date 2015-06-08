%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bounds1D.Discretize}
\alias{bounds1D.Discretize}
\alias{Discretize.bounds1D}
\alias{bounds1D,Discretize-method}

\title{Min, mean and max density on a 1D cube}

\description{
  Internal function.
}

\usage{
\method{bounds1D}{Discretize}(this, cube, mu, sigma2, len=100, ...)
}

\arguments{
 \item{cube}{ The cube under consideration which is a (2x2) matrix containing the bounds of the variables (columnwise) (bivariate case).   }
 \item{mu}{ The mean. }
 \item{sigma2}{ The variance. }
 \item{len}{ The number of samples of each coordinate.}
 \item{...}{Not used.}
}

\value{
  Min, mean and max density on a 1D cube.
}

\author{Lars Relund \email{lars@relund.dk}}

\seealso{
  For more information see \code{\link{Discretize}}.
}


\keyword{internal}
\keyword{methods}
