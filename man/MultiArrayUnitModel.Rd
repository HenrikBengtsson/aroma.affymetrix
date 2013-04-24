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
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{MultiArrayUnitModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineCnPlm}}, \emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{AvgCnPlm}}, \emph{\link[aroma.affymetrix]{AvgPlm}}, \emph{\link[aroma.affymetrix]{AvgSnpPlm}}, \emph{\link[aroma.affymetrix]{ExonRmaPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddCnPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public abstract static class \bold{MultiArrayUnitModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



 This abstract class represents a unit model that fits one model per unit
 based on signals for all arrays in the data set.
 The nature of a multi-array unit model is that all arrays must be
 available at the time of the fit and the estimated parameters will
 depend on the data from all arrays.
 Thus, if the signals in one array changes the model has to be refitted.
}

\usage{
MultiArrayUnitModel(..., listOfPriors=NULL)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{UnitModel}}.}
  \item{listOfPriors}{A \code{\link[base]{list}} of priors to be used when fitting
   the model.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getListOfPriors} \tab  -\cr
 \tab \code{setListOfPriors} \tab  -\cr
}


 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction, getParameters

 \bold{Methods inherited from Model}:\cr
as.character, fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getParameterSet, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
