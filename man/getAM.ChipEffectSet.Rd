%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.xam.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getAM.ChipEffectSet}
\alias{getAM.ChipEffectSet}
\alias{ChipEffectSet.getAM}
\alias{getAM.ChipEffectSet}
\alias{getAM,ChipEffectSet-method}


\title{Gets the log-intensities and log-ratios of chip effects of the set relative to a reference chip effect file}

\description{
 Gets the log-intensities and log-ratios of chip effects of the set relative to a reference chip effect file.
}

\usage{\method{getAM}{ChipEffectSet}(this, other, units=NULL, ..., verbose=FALSE)}

\arguments{
  \item{other}{The second \code{\link{ChipEffectFile}} object used as the
    reference.}
  \item{units}{(The subset of units to be matched.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an Jx2xK \code{\link[base]{array}} where J is the number of units, and K is
 the number of arrays (arrays are always the last dimension).
 The names of the columns are A (log-intensities) and M (log-ratios).
 The names of the rows are the unit indices (as indexed by the CDF).
 The rows are ordered according to \code{units} arguments.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
