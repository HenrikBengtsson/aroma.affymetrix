%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.xam.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getXAM.ChipEffectSet}
\alias{getXAM.ChipEffectSet}
\alias{ChipEffectSet.getXAM}
\alias{getXAM,ChipEffectSet-method}

\title{Gets the physical position, log-intensities and log-ratios of chip effects for two arrays}

\description{
 Gets the physical position, log-intensities and log-ratios of chip effects for two arrays of units on a certain chromosome.
}

\usage{
\method{getXAM}{ChipEffectSet}(this, other, chromosome, units=NULL, ..., verbose=FALSE)
}

\arguments{
  \item{other}{The second \code{\link{ChipEffectFile}} object used as the
    reference.}
  \item{chromosome}{(The chromosome for which results should be returned.}
  \item{units}{(The subset of units to be matched.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a Nx3 matrix where N is the number of units returned.
 The names of the columns are X (physical position in a given chromosome),
 A (log-intensities) and M (log-ratios).
 The names of the rows are the unit indices (as indexed by the CDF).
 \emph{Note: The rows are ordered according to chromosomal position.}
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:getAM.ChipEffectSet]{*getAM}()}.
  For more information see \code{\link{ChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
