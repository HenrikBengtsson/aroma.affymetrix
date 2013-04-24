%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{isDuplicated.AffymetrixCelSet}
\alias{isDuplicated.AffymetrixCelSet}
\alias{AffymetrixCelSet.isDuplicated}
\alias{isDuplicated,AffymetrixCelSet-method}

\title{Identifies duplicated CEL files}

\description{
  Identifies duplicated CEL files by comparing the timestamps in the CEL headers.
}

\usage{
\method{isDuplicated}{AffymetrixCelSet}(this, ..., verbose=FALSE)
}

\arguments{
 \item{...}{Not used.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns a \code{\link[base]{logical}} \code{\link[base]{vector}} of length equal to the number of files
  in the set.
  An element with value \code{\link[base:logical]{TRUE}} indicates that the corresponding CEL file
  has the same time stamp as another preceeding CEL file.
}

\author{Henrik Bengtsson}

\seealso{
  Internally \code{\link[base]{duplicated}}() is used to compare timestamps.
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
