%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnagCfhFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnits.CnagCfhFile}
\alias{readUnits.CnagCfhFile}
\alias{CnagCfhFile.readUnits}
\alias{readUnits.CnagCfhFile}
\alias{readUnits,CnagCfhFile-method}


\title{Reads CEL data unit by unit}

\description{
 Reads CEL data unit by unit for all or a subset of units (probeset).
}

\usage{\method{readUnits}{CnagCfhFile}(this, units=NULL, ..., verbose=FALSE)}

\arguments{
  \item{units}{The units to be read. If \code{\link[base]{NULL}}, all units are read.}
  \item{cdf}{An alternative CDF structure to be used.  This overrides
    the \code{units} arguments.}
  \item{...}{Arguments passed to \code{getUnits()} of the
    \code{\link{AffymetrixCdfFile}} class (if \code{cdf} was not specified),
    but also to the \code{\link[affxparser]{readCelUnits}} methods.}
}

\value{
 Returns the \code{\link[base]{list}} structure that \code{\link[affxparser]{readCelUnits}} returns.
}

\section{Caching}{
  CEL data is neither cached in memory nor on file by this method.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:updateUnits.CnagCfhFile]{*updateUnits}()}.
  For more information see \code{\link{CnagCfhFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
