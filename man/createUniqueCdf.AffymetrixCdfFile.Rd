%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.UNIQUE.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{createUniqueCdf.AffymetrixCdfFile}
\alias{createUniqueCdf.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.createUniqueCdf}
\alias{createUniqueCdf.AffymetrixCdfFile}
\alias{createUniqueCdf,AffymetrixCdfFile-method}

\alias{AffymetrixCdfFile.createUnique}
\alias{createUnique.AffymetrixCdfFile}
\alias{createUnique,AffymetrixCdfFile-method}


\title{Creates a unique-cell version of the CDF}

\description{
 In some cases, single probes map to multiple genome locations.  In cases
 where you may later want to store a probe estimate (e.g. a probe effect
 or a residual), you will not be able to store more than one per cell.
 The unique CDF facilitates this by making the cell indices unique
 (essentially copying the multimapping probes).
}

\usage{\method{createUniqueCdf}{AffymetrixCdfFile}(this, chipType=getChipType(this), tags="unique", sep=",", path=NULL, units=NULL, ..., ram=NULL, verbose=TRUE)}

\arguments{
  \item{chipType}{The chip type of the new CDF.}
  \item{tags}{Tags added to the chip type of the new CDF.}
  \item{sep}{A string separating the chip type and the tags string.}
  \item{path}{The path where to store the new CDF file.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{writeCdf}}.}
  \item{ram}{A \code{\link[base]{double}} saying if more or less units should be converted
     per chunk.  A smaller value uses less memory.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the unique CDF as an \code{\link{AffymetrixCdfFile}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
