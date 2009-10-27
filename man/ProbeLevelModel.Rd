%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbeLevelModel}
\docType{class}
\alias{ProbeLevelModel}


 \title{The ProbeLevelModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbeLevelModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MultiArrayUnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{ProbeLevelModel}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffineCnPlm}, \link[aroma.affymetrix]{AffinePlm}, \link[aroma.affymetrix]{AffineSnpPlm}, \link[aroma.affymetrix]{AvgCnPlm}, \link[aroma.affymetrix]{AvgPlm}, \link[aroma.affymetrix]{AvgSnpPlm}, \link{ExonRmaPlm}, \link[aroma.affymetrix]{HetLogAddCnPlm}, \link[aroma.affymetrix]{HetLogAddPlm}, \link[aroma.affymetrix]{HetLogAddSnpPlm}, \link[aroma.affymetrix]{MbeiCnPlm}, \link[aroma.affymetrix]{MbeiPlm}, \link[aroma.affymetrix]{MbeiSnpPlm}, \link[aroma.affymetrix]{RmaCnPlm}, \link{RmaPlm}, \link[aroma.affymetrix]{RmaSnpPlm}\cr

public static class \bold{ProbeLevelModel}\cr
extends \emph{\link[aroma.affymetrix]{MultiArrayUnitModel}}\cr



  This abstract class represents a probe-level model (PLM) as defined
  by the \pkg{affyPLM} package:
    "A [...] PLM is a model that is fit to probe-intensity data.
     More specifically, it is where we fit a model with probe level
     and chip level parameters on a probeset by probeset basis",
  where the more general case for a probeset is a \emph{unit group}
  in Affymetrix CDF terms.
 }

 \usage{ProbeLevelModel(..., standardize=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{MultiArrayUnitModel}}.}
   \item{standardize}{If \code{\link[base:logical]{TRUE}}, chip-effect and probe-affinity estimates are
      rescaled such that the product of the probe affinities is one.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{calculateResidualSet} \tab  -\cr
 \tab \code{calculateWeights} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.ProbeLevelModel]{fit}} \tab Estimates the model parameters.\cr
 \tab \code{\link[aroma.affymetrix:getChipEffectSet.ProbeLevelModel]{getChipEffectSet}} \tab Gets the set of chip effects for this model.\cr
 \tab \code{\link[aroma.affymetrix:getProbeAffinityFile.ProbeLevelModel]{getProbeAffinityFile}} \tab Gets the probe affinities for this model.\cr
 \tab \code{getResidualSet} \tab  -\cr
 \tab \code{getWeightsSet} \tab  -\cr
}


 \bold{Methods inherited from MultiArrayUnitModel}:\cr
getListOfPriors, setListOfPriors, validate

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \details{
   In order to minimize the risk for mistakes, but also to be able compare
   results from different PLMs, all PLM subclasses must meet the following
   criteria:
   \enumerate{
     \item All parameter estimates must be (stored and returned) on the
       intensity scale, e.g. log-additive models such as \code{\link{RmaPlm}}
       have to transform the parameters on the log-scale to the intensity
       scale.
     \item The probe-affinity estimates \eqn{\phi_k} for a unit group
       must be constrained such that \eqn{\prod_k \phi_k = 1},
       or equivalently if \eqn{\phi_k > 0},\eqn{\sum_k \log(\phi_k) = 0}.
   }
   Note that the above probe-affinity constraint guarantees that the
   estimated chip effects across models are on the same scale.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more details on probe-level models, please see
   the \pkg{affyPLM} package.
 }
\keyword{classes}
