%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RmaPlm}
\docType{class}
\alias{RmaPlm}


 \title{The RmaPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RmaPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MultiArrayUnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{RmaPlm}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ExonRmaPlm}, \link[aroma.affymetrix]{HetLogAddCnPlm}, \link[aroma.affymetrix]{HetLogAddPlm}, \link[aroma.affymetrix]{HetLogAddSnpPlm}, \link[aroma.affymetrix]{RmaCnPlm}, \link[aroma.affymetrix]{RmaSnpPlm}\cr

public static class \bold{RmaPlm}\cr
extends \link[aroma.affymetrix]{ProbeLevelModel}\cr



  This class represents the log-additive model part of the Robust Multichip
  Analysis (RMA) method described in Irizarry et al (2003).
 }

 \usage{RmaPlm(..., flavor=c("affyPLM", "oligo"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
   \item{flavor}{A \code{\link[base]{character}} string specifying what model fitting algorithm
     to be used.  This makes it possible to get identical estimates as other
     packages.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getAsteriskTags} \tab  -\cr
}


 \bold{Methods inherited from ProbeLevelModel}:\cr
calculateResidualSet, calculateWeights, fit, getAsteriskTags, getCalculateResidualsFunction, getChipEffectSet, getProbeAffinityFile, getResidualSet, getWeightsSet

 \bold{Methods inherited from MultiArrayUnitModel}:\cr
getListOfPriors, setListOfPriors, validate

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Model}{
   For a single unit group, the log-additive model of RMA is:

    \deqn{log_2(y_{ik}) = \beta_i + \alpha_k + \varepsilon_{ik}}

   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I},
   and \eqn{\alpha_k} are the probe affinities for probes \eqn{k=1,...,K}.
   The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
   The model is constrained such that \eqn{\sum_k{\alpha_k} = 0}.

   Note that all PLM classes must return parameters on the intensity scale.
   For this class that means that \eqn{\theta_i = 2^\beta_i} and
   \eqn{\phi_k = 2^\alpha_k} are returned.
 }

 \section{Different flavors of model fitting}{
   There are a few differ algorithms available for fitting the same
   probe-level model.  The default and recommended method
   (\code{flavor="affyPLM"}) uses the implementation in the
   \pkg{preprocessCore} package which fits the model parameters robustly
   using an M-estimator (the method used to be in \pkg{affyPLM}).

   Alternatively, other model-fitting algorithms are available.
   The algorithm (\code{flavor="oligo"}) used by the \pkg{oligo} package,
   which originates from the \pkg{affy} packages, fits the model using
   median polish, which is a non-robust estimator.  Note that this algorithm
   does not constraint the probe-effect parameters to multiply to one on
   the intensity scale.  Since the internal function does not return these
   estimates, we can neither rescale them.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}.
  NAR, 2003, 31, e15.\cr
 }
\keyword{classes}
