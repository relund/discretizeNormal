%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bounds2D.Discretize}
\alias{bounds2D.Discretize}
\alias{Discretize.bounds2D}
\alias{bounds2D,Discretize-method}

\title{Min, mean and max density on a 2D cube}

\description{
  Internal function.
}

\usage{
\method{bounds2D}{Discretize}(this, cube, mu, sigma, len=100, ...)
}

\arguments{
 \item{cube}{ The cube under consideration which is a (2x2) matrix containing the bounds of the variables (columnwise) (bivariate case).   }
 \item{mu}{ The mean. }
 \item{sigma}{ The covariate matrix. }
 \item{len}{ The number of samples of each coordinate.}
 \item{...}{Not used.}
}

\value{
  Min, mean and max density on a 2D cube.
}

\author{Lars Relund \email{lars@relund.dk}}

\seealso{
  For more information see \code{\link{Discretize}}.
}


\keyword{internal}
\keyword{methods}
