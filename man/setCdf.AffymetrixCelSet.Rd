%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setCdf.AffymetrixCelSet}
\alias{setCdf.AffymetrixCelSet}
\alias{AffymetrixCelSet.setCdf}
\alias{setCdf,AffymetrixCelSet-method}

\title{Sets the CDF structure for this CEL set}

\description{
 Sets the CDF structure for this CEL set.  This structures is assigned to all CEL files in the set.
}

\usage{
\method{setCdf}{AffymetrixCelSet}(this, cdf, verbose=FALSE, ..., .checkArgs=TRUE)
}

\arguments{
  \item{cdf}{An \code{\link{AffymetrixCdfFile}} object.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
    May also be a \code{\link[R.utils]{Verbose}} object.}
  \item{...}{Not used.}
  \item{.checkArgs}{(Internal) If \code{\link[base:logical]{FALSE}}, arguments are not validated.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:getCdf.AffymetrixCelSet]{*getCdf}()}.
  For more information see \code{\link{AffymetrixCelSet}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
