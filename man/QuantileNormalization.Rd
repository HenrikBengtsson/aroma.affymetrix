%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantileNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QuantileNormalization}
\docType{class}
\alias{QuantileNormalization}


\title{The QuantileNormalization class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class QuantileNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{QuantileNormalization}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{DChipQuantileNormalization}\cr

public static class \bold{QuantileNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



 This class represents a normalization function that transforms the
 probe-level signals towards the same empirical distribution.
}

\usage{
QuantileNormalization(..., subsetToUpdate=NULL, typesToUpdate=NULL,
  targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate)
}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{ProbeLevelTransform}}.}
  \item{subsetToUpdate}{The probes to be updated.
    If \code{\link[base]{NULL}}, all probes are updated.}
  \item{typesToUpdate}{Types of probes to be updated.}
  \item{targetDistribution}{A \code{\link[base]{numeric}} \code{\link[base]{vector}}.  The empirical
    distribution to which all arrays should be normalized to.}
  \item{subsetToAvg}{The probes to calculate average empirical
    distribution over.  If a single \code{\link[base]{numeric}} in (0,1), then this
    fraction of all probes will be used.
    If \code{\link[base]{NULL}}, all probes are considered.}
  \item{typesToAvg}{Types of probes to be used when calculating the
    average empirical distribution.
    If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
    probes are used, respectively. If \code{"pmmm"} both types are used.
  }
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.QuantileNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from ProbeLevelTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, asThis


}

\examples{\dontrun{
  for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find any dataset
path <- NULL
if (is.null(path))
  break

ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
norm <- QuantileNormalization(ds, subsetToAvg=1/3)
dsQN <- process(norm, verbose=verbose)
print(dsQN)

} # for (zzz in 0)
rm(zzz)

}}

\author{Henrik Bengtsson}
\keyword{classes}
