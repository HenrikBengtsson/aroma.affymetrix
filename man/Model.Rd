%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Model.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Model}
\docType{class}
\alias{Model}


\title{The Model class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class Model}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{Model}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineCnPlm}}, \emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{AlleleSummation}}, \emph{\link[aroma.affymetrix]{AvgCnPlm}}, \emph{\link[aroma.affymetrix]{AvgPlm}}, \emph{\link[aroma.affymetrix]{AvgSnpPlm}}, \link[aroma.affymetrix]{CrlmmModel}, \emph{\link[aroma.affymetrix]{ExonRmaPlm}}, \emph{\link[aroma.affymetrix]{FirmaModel}}, \emph{\link[aroma.affymetrix]{HetLogAddCnPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{MultiArrayUnitModel}}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}, \emph{\link[aroma.affymetrix]{SingleArrayUnitModel}}, \emph{\link[aroma.affymetrix]{UnitModel}}\cr

public abstract static class \bold{Model}\cr
extends \link[aroma.core]{ParametersInterface}\cr



 This class is abstract and represents a generic model that applies
 to a data set.
}

\usage{
Model(dataSet=NULL, tags="*", ..., .onUnknownArgs=c("error", "warning", "ignore"))
}

\arguments{
  \item{dataSet}{The data set to which this model should be fitted.}
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
     the input data set.}
  \item{...}{Not used.}
  \item{.onUnknownArgs}{A \code{\link[base]{character}} string specifying what should occur
     if there are unknown arguments in \code{...}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fit} \tab  -\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
