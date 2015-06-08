%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{klBound2D.Discretize}
\alias{klBound2D.Discretize}
\alias{Discretize.klBound2D}
\alias{klBound2D,Discretize-method}

\title{Upper bound on KL distance on a 2D cube}

\description{
  Internal function.
}

\usage{
\method{klBound2D}{Discretize}(this, cube, mu, sigma, ...)
}

\arguments{
 \item{cube}{ The cube under consideration which is a (2x2) matrix containing the bounds of the variables (columnwise) (bivariate case).   }
 \item{mu}{ The mean. }
 \item{sigma}{ The covariate matrix. }
 \item{...}{Not used.}
}

\value{
  Upper bound on KL distance on a 2D cube.
}

\author{Lars Relund \email{lars@relund.dk}}

\seealso{
  For more information see \code{\link{Discretize}}.
}


\keyword{internal}
\keyword{methods}
