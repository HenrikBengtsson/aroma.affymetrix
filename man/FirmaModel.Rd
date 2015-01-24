%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FirmaModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FirmaModel}
\docType{class}
\alias{FirmaModel}


\title{The FirmaModel class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class FirmaModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{FirmaModel}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{FirmaModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



 This class represents the FIRMA (Finding Isoforms using RMA) alternative
 splicing model.

}

\usage{
FirmaModel(rmaPlm=NULL, summaryMethod=c("median", "upperQuartile", "max"),
  operateOn=c("residuals", "weights"), ...)
}

\arguments{
  \item{rmaPlm}{An @RmaPlm object.}
  \item{summaryMethod}{A \code{\link[base]{character}} specifying what summarization method should be used.}
  \item{operateOn}{A \code{\link[base]{character}} specifying what statistic to operate on.}
  \item{...}{Arguments passed to constructor of \code{\link{UnitModel}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fit} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{getFirmaSet} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPlm} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
}


 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction, getParameters

 \bold{Methods inherited from Model}:\cr
as.character, fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Ken Simpson, Henrik Bengtsson}
\keyword{classes}
