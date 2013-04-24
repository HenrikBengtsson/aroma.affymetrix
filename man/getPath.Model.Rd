%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Model.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getPath.Model}
\alias{getPath.Model}
\alias{Model.getPath}
\alias{getPath,Model-method}

\title{Gets the path of this model}

\description{
 Gets the path of this model where the parameter files are located.
 The directory is created, if missing.
}

\usage{
\method{getPath}{Model}(this, ...)
}

\arguments{
 \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{character}} string.
}

\details{
 If the path does not exist, it is created.
}

\seealso{
  For more information see \code{\link{Model}}.
}
\keyword{internal}
\keyword{methods}
