%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.SNPs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{nbrOfSnps.AffymetrixCdfFile}
\alias{nbrOfSnps.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.nbrOfSnps}
\alias{nbrOfSnps,AffymetrixCdfFile-method}

\title{Gets the number of SNPs}

\description{
  Gets the number of SNPs.
}

\usage{
\method{nbrOfSnps}{AffymetrixCdfFile}(this, ...)
}

\arguments{
 \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:getSnpNames.AffymetrixCdfFile]{*getSnpNames}()}.}
}

\value{
  Returns an \code{\link[base]{integer}}.
}

\author{Henrik Bengtsson}

\seealso{
  Internally, \code{\link[aroma.affymetrix:getSnpNames.AffymetrixCdfFile]{*getSnpNames}()} is used to identify SNPs.

  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
