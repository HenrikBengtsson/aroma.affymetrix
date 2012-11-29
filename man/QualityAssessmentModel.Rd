%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QualityAssessmentModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QualityAssessmentModel}
\docType{class}
\alias{QualityAssessmentModel}


\title{The QualityAssessmentModel class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class QualityAssessmentModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{QualityAssessmentModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{QualityAssessmentModel}\cr
extends \link[R.oo]{Object}\cr


}

\usage{QualityAssessmentModel(plm=NULL, tags="*", ...)}

\arguments{
  \item{plm}{A \code{\link{ProbeLevelModel}}.}
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getChipEffectSet} \tab  -\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getPlm} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getResiduals.QualityAssessmentModel]{getResiduals}} \tab Calculates the residuals from a probe-level model.\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getWeights.QualityAssessmentModel]{getWeights}} \tab Calculates the weights from the robust fit to a probe-level model.\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{plotNuse} \tab  -\cr
 \tab \code{plotRle} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}
\keyword{classes}
