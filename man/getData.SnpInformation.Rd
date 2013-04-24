%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getData.SnpInformation}
\alias{getData.SnpInformation}
\alias{SnpInformation.getData}
\alias{getData,SnpInformation-method}

\title{Gets all or a subset of the genome information data}

\description{
  Gets all or a subset of the genome information data.
}

\usage{
\method{getData}{SnpInformation}(this, units=NULL, fields=c("fragmentLength", "start", "stop"), orderBy=NULL, ..., force=FALSE, verbose=FALSE)
}

\arguments{
 \item{units}{The units for which the data should be returned.}
 \item{fields}{The fields to be returned.}
 \item{orderBy}{The fields by which the returned data frame should be
     ordered.}
 \item{...}{Named arguments used to select a subset of the units to be
     returned.  Either a value to be compared to or a \code{\link[base]{function}} returning
     \code{\link[base:logical]{TRUE}} or \code{\link[base:logical]{FALSE}}.}
}

\value{
  Returns a \code{\link[base]{data.frame}}, where the row names correspond to unit indices
  as defined by the CDF.
}

\seealso{
  For more information see \code{\link{SnpInformation}}.
}
\keyword{internal}
\keyword{methods}
