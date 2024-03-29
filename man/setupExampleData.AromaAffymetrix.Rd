%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  setupExampleData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{setupExampleData.AromaAffymetrix}
\alias{setupExampleData.AromaAffymetrix}
\alias{AromaAffymetrix.setupExampleData}
\alias{setupExampleData,AromaAffymetrix-method}
\alias{setupExampleData}

\title{Setups example data in the current directory}

\description{
 Setups example data in the current directory.
}

\usage{
\method{setupExampleData}{AromaAffymetrix}(pkg, dataset=NULL, chipType=NULL, dirs=c("annotationData", "rawData"),
  mustWork=TRUE, validate=FALSE, ...)
}

\arguments{
  \item{dataset, chipType}{\code{\link[base]{character}} strings specifying the data set
   and the chip type to install.}
  \item{dirs}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying which directories to setup.}
  \item{mustWork}{If \code{\link[base:logical]{TRUE}}, an error is thrown if the requested data set
   could not be installed, otherwise not.}
  \item{validate}{If \code{\link[base:logical]{TRUE}}, the installed files are also validated,
   otherwise not.}
  \item{...}{Not used.}
}

\value{
  Returns (invisibly) \code{\link[base:logical]{TRUE}} if all requested data was installed,
  otherwise \code{\link[base:logical]{FALSE}}.
}

\author{Henrik Bengtsson}


\keyword{internal}
\keyword{methods}
