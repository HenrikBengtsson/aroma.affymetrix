%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  randomSeed.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{randomSeed}
\alias{randomSeed}


\title{Sets and resets the .Random.seed in the global environment}

\description{
 Sets and resets the .Random.seed in the global environment.
}

\usage{
 randomSeed(action=c("set", "advance", "reset", "get"), seed=NULL, kind=NULL, backup=TRUE)
}

\arguments{
  \item{action}{A \code{\link[base]{character}} string specifying the action.}
  \item{seed}{Random seed to be set; only for \code{action="set"}.
    If \code{length(seed) == 1}, then \code{set.seed(seed)} is
    used, otherwise \code{.Random.seed} is assigned the value.}
  \item{kind}{(optional) A \code{\link[base]{character}} string specifying type of
    random number generator to use, cf. \code{\link[base]{RNGkind}}().}
  \item{backup}{If \code{\link[base:logical]{TRUE}}, the previous (seed, kind) state is recorded
    such that it can be reset later.}
}

\value{
  Returns a \code{.Random.seed}.
}

\author{Henrik Bengtsson}


\keyword{internal}
