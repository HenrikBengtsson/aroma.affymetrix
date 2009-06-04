%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.calculateBaseline.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calculateBaseline.ChipEffectSet}
\alias{calculateBaseline.ChipEffectSet}
\alias{ChipEffectSet.calculateBaseline}
\alias{calculateBaseline.ChipEffectSet}
\alias{calculateBaseline,ChipEffectSet-method}


\title{Estimates the baseline signal chromosome by chromosome}

\description{
 Estimates the baseline signal chromosome by chromosome.
}

\usage{\method{calculateBaseline}{ChipEffectSet}(this, chromosomes=NULL, ploidy=2, defaultPloidy=NA, all=FALSE, force=FALSE, verbose=FALSE, ...)}

\arguments{
  \item{chromosomes}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying for which chromsosomes
    the baseline should be estimated.
    If \code{\link[base]{NULL}}, all chromosomes are considered.}
  \item{ploidy}{An \code{\link[base]{integer}} specifying the ploidy that the baseline
    should have.}
  \item{defaultPloidy}{An \code{\link[base]{integer}} specifying the default ploidy of
    chromosomes that have not explicitly been allocated one.}
  \item{all}{If \code{\link[base:logical]{TRUE}}, signals are averaged also for cells that are not
    on the genome.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, the CEL file that stores the is recreated.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link{getAverageFile}}().
  For more information see \code{\link{ChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
