%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getChipEffectSet.ProbeLevelModel}
\alias{getChipEffectSet.ProbeLevelModel}
\alias{ProbeLevelModel.getChipEffectSet}
\alias{getChipEffectSet,ProbeLevelModel-method}
\alias{ProbeLevelModel.getChipEffects}
\alias{getChipEffects.ProbeLevelModel}
\alias{getChipEffects,ProbeLevelModel-method}


\title{Gets the set of chip effects for this model}

\description{
 Gets the set of chip effects for this model.
 There is one chip-effect file per array.
}

\usage{
\method{getChipEffectSet}{ProbeLevelModel}(this, ..., verbose=FALSE)
}

\arguments{
  \item{...}{Arguments passed to \code{getMonocellCdf()} of
                                                 \code{\link{AffymetrixCdfFile}}.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link{ChipEffectSet}} object.
}

\seealso{
  For more information see \code{\link{ProbeLevelModel}}.
}
\keyword{internal}
\keyword{methods}
