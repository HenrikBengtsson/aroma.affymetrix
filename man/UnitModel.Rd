%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{UnitModel}
\docType{class}
\alias{UnitModel}


\title{The UnitModel class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class UnitModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{UnitModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineCnPlm}}, \emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{AlleleSummation}}, \emph{\link[aroma.affymetrix]{AvgCnPlm}}, \emph{\link[aroma.affymetrix]{AvgPlm}}, \emph{\link[aroma.affymetrix]{AvgSnpPlm}}, \emph{\link[aroma.affymetrix]{ExonRmaPlm}}, \emph{\link[aroma.affymetrix]{FirmaModel}}, \emph{\link[aroma.affymetrix]{HetLogAddCnPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{MultiArrayUnitModel}}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}, \emph{\link[aroma.affymetrix]{SingleArrayUnitModel}}\cr

public abstract static class \bold{UnitModel}\cr
extends \emph{\link[aroma.affymetrix]{Model}}\cr



 This class is abstract and represents a generic unit model, i.e.
 a model that is applied to each unit separately.
}

\usage{
UnitModel(dataSet=NULL, probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)", "pm+mm"),
  shift=0, ...)
}

\arguments{
  \item{dataSet}{An \code{\link{AffymetrixCelSet}} to which this model should
     be fitted.}
  \item{probeModel}{A \code{\link[base]{character}} string specifying how PM and MM values
     should be modelled.  By default only PM signals are used.}
  \item{shift}{An optional amount the signals should be shifted
     (translated) before fitting the model.}
  \item{...}{Arguments passed to the constructor of \code{\link{Model}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
}


 \bold{Methods inherited from Model}:\cr
as.character, fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
