%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AvgPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AvgPlm}
\docType{class}
\alias{AvgPlm}


\title{The AvgPlm class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AvgPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MultiArrayUnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AvgPlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AvgCnPlm}}, \emph{\link[aroma.affymetrix]{AvgSnpPlm}}\cr

public abstract static class \bold{AvgPlm}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



 This class represents a PLM where the probe intensities are averaged
 assuming identical probe affinities.
 For instance, one may assume that replicated probes with identical
 sequences have the same probe affinities, cf. the GenomeWideSNP\_6
 chip type.
}

\usage{
AvgPlm(..., flavor=c("median", "mean"))
}

\arguments{
  \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
  \item{flavor}{A \code{\link[base]{character}} string specifying what model fitting algorithm
    to be used.  This makes it possible to get identical estimates as other
    packages.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from ProbeLevelModel}:\cr
calculateResidualSet, calculateWeights, fit, getAsteriskTags, getCalculateResidualsFunction, getChipEffectSet, getProbeAffinityFile, getResidualSet, getRootPath, getWeightsSet

 \bold{Methods inherited from MultiArrayUnitModel}:\cr
getListOfPriors, setListOfPriors, validate

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction, getParameters

 \bold{Methods inherited from Model}:\cr
as.character, fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getParameterSet, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, asThis


}

\section{Model}{
  For a single unit group, the averaging PLM of K probes is:

   \deqn{y_{ik} = \theta_i + \varepsilon_{ik}}

  where \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I}.
  The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
}

\section{Different flavors of model fitting}{
  The above model can be fitted in two ways, either robustly or
  non-robustly.
  Use argument \code{flavor="mean"} to fit the model non-robustly, i.e.

   \deqn{\hat{\theta}_{i} = 1/K \sum_k y_{ik}}.

  Use argument \code{flavor="median"} to fit the model robustly, i.e.

   \deqn{\hat{\theta}_{i} = median_k y_{ik}}.

  Missing values are always excluded.
}

\author{Henrik Bengtsson}
\keyword{classes}
