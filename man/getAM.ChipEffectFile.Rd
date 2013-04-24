%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.xam.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getAM.ChipEffectFile}
\alias{getAM.ChipEffectFile}
\alias{ChipEffectFile.getAM}
\alias{getAM,ChipEffectFile-method}

\title{Gets the log-intensities and log-ratios of chip effects for two arrays}

\description{
 Gets the log-intensities and log-ratios of chip effects for two arrays.
}

\usage{
\method{getAM}{ChipEffectFile}(this, other, units=NULL, ..., verbose=FALSE)
}

\arguments{
  \item{other}{The second \code{\link{ChipEffectFile}} object used as the
    reference.}
  \item{units}{(The subset of units to be matched.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a Nx2 matrix where N is the number of units returned.
 The names of the columns are A (log-intensities) and M (log-ratios).
 The names of the rows are the unit indices (as indexed by the CDF).
 The rows are ordered according to \code{units} arguments.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:getXAM.ChipEffectFile]{*getXAM}()}.
  For more information see \code{\link{ChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
