%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnagCfhSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{isDuplicated.CnagCfhSet}
\alias{isDuplicated.CnagCfhSet}
\alias{CnagCfhSet.isDuplicated}
\alias{isDuplicated,CnagCfhSet-method}

\title{Identifies duplicated CFH files}

\description{
  Identifies duplicated CFH files by comparing the timestamps in the CFH headers.
}

\usage{
\method{isDuplicated}{CnagCfhSet}(this, ...)
}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{logical}} \code{\link[base]{vector}} of length equal to the number of files
  in the set.
  An element with value \code{\link[base:logical]{TRUE}} indicates that the corresponding CFH file
  has the same time stamp as another preceeding CFH file.
}

\seealso{
  Internally \code{\link[base]{duplicated}}() is used to compare timestamps.
  For more information see \code{\link{CnagCfhSet}}.
}
\keyword{internal}
\keyword{methods}
