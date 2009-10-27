%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiArrayUnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MultiArrayUnitModel}
\docType{class}
\alias{MultiArrayUnitModel}


 \title{The MultiArrayUnitModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class MultiArrayUnitModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{MultiArrayUnitModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffineCnPlm}, \link[aroma.affymetrix]{AffinePlm}, \link[aroma.affymetrix]{AffineSnpPlm}, \link[aroma.affymetrix]{AvgCnPlm}, \link[aroma.affymetrix]{AvgPlm}, \link[aroma.affymetrix]{AvgSnpPlm}, \link{ExonRmaPlm}, \link[aroma.affymetrix]{HetLogAddCnPlm}, \link[aroma.affymetrix]{HetLogAddPlm}, \link[aroma.affymetrix]{HetLogAddSnpPlm}, \link[aroma.affymetrix]{MbeiCnPlm}, \link[aroma.affymetrix]{MbeiPlm}, \link[aroma.affymetrix]{MbeiSnpPlm}, \link[aroma.affymetrix]{ProbeLevelModel}, \link[aroma.affymetrix]{RmaCnPlm}, \link{RmaPlm}, \link[aroma.affymetrix]{RmaSnpPlm}\cr

public abstract static class \bold{MultiArrayUnitModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



  This abstract class represents a unit model that fits one model per unit
  based on signals for all arrays in the data set.
  The nature of a multi-array unit model is that all arrays must be
  available at the time of the fit and the estimated parameters will
  depend on the data from all arrays.
  Thus, if the signals in one array changes the model has to be refitted.
 }

 \usage{MultiArrayUnitModel(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{UnitModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getListOfPriors} \tab  -\cr
 \tab \code{setListOfPriors} \tab  -\cr
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
