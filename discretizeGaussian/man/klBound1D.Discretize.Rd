%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  discretize.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{klBound1D.Discretize}
\alias{klBound1D.Discretize}
\alias{Discretize.klBound1D}
\alias{klBound1D,Discretize-method}

\title{Upper bound on KL distance on a 1D cube}

\description{
  Internal function.
}

\usage{
\method{klBound1D}{Discretize}(this, cube, mu, sigma2, ...)
}

\arguments{
 \item{cube}{ The cube under consideration which is vector of length 2 containing the upper and lower bounds (univariate case). }
 \item{mu}{ The mean. }
 \item{sigma2}{ The variance. }
 \item{...}{ Not used. }
}

\value{
  Upper bound on KL distance on a 1D cube.
}

\author{Lars Relund \email{lars@relund.dk}}

\seealso{
  For more information see \code{\link{Discretize}}.
}


\keyword{internal}
\keyword{methods}
