%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FragmentLengthNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FragmentLengthNormalization}
\docType{class}
\alias{FragmentLengthNormalization}


\title{The FragmentLengthNormalization class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class FragmentLengthNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectTransform}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{FragmentLengthNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FragmentLengthNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectTransform}}\cr



 This class represents a normalization method that corrects for PCR
 fragment length effects on copy-number chip-effect estimates.
}

\usage{
FragmentLengthNormalization(dataSet=NULL, ..., target=targetFunctions, subsetToFit="-XY",
  lengthRange=NULL, onMissing=c("median", "ignore"), shift=0, targetFunctions=NULL)
}

\arguments{
  \item{dataSet}{A \code{\link{SnpChipEffectSet}}.}
  \item{...}{Additional arguments passed to the constructor of
    \code{\link{ChipEffectTransform}}.}
  \item{target}{(Optional) A \code{\link[base]{character}} string or a list of \code{\link[base]{function}}s
    specifying what to normalize toward.
    For each enzyme there is one target function to which all arrays
    should be normalized to.}
  \item{subsetToFit}{The units from which the normalization curve should
    be estimated.  If \code{\link[base]{NULL}}, all are considered.}
  \item{lengthRange}{If given, a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length 2 specifying
    the range of fragment lengths considered.  All fragments with lengths
    outside this range are treated as if they were missing.}
  \item{onMissing}{Specifies how to normalize units for which the
    fragment lengths are unknown.}
  \item{shift}{An optional amount the data points should be shifted
     (translated).}
  \item{targetFunctions}{Deprecated.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCdf} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, findFilesTodo, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\section{Requirements}{
  This class requires a SNP information annotation file for the
  chip type to be normalized.
}

\details{
  For SNPs, the normalization function is estimated based on the total
  chip effects, i.e. the sum of the allele signals.  The normalizing
  is done by rescale the chip effects on the intensity scale such that
  the mean of the total chip effects are the same across samples for
  any given fragment length.  For allele-specific estimates, both alleles
  are always rescaled by the same amount.  Thus, when normalizing
  allele-specific chip effects, the total chip effects is change, but not
  the relative allele signal, e.g. the allele B frequency.
}

\author{Henrik Bengtsson}
\keyword{classes}
