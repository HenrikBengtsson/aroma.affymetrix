%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AlleleSummation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AlleleSummation}
\docType{class}
\alias{AlleleSummation}


 \title{The AlleleSummation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AlleleSummation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{AlleleSummation}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AlleleSummation}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



  This class takes allele-specific chip effect estimates of a
  SnpChipEffectSet and returns a CnChipEffectSet holding the summed
  allele estimates.
 }

 \usage{AlleleSummation(dataSet=NULL, ignoreNAs=TRUE, ...)}

 \arguments{
   \item{dataSet}{A \code{\link{SnpChipEffectSet}}.}
   \item{ignoreNAs}{If \code{\link[base:logical]{TRUE}}, missing values are excluded when summing
      the signals from the two alleles.}
   \item{...}{Arguments passed to \code{\link{UnitModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipEffectSet.AlleleSummation]{getChipEffectSet}} \tab Gets the set of chip effects for this model.\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
