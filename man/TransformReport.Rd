%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TransformReport.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{TransformReport}
\docType{class}
\alias{TransformReport}


\title{The TransformReport class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class TransformReport}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{TransformReport}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TransformReport}\cr
extends \link[R.oo]{Object}\cr


}

\usage{
TransformReport(inSet=NULL, outSet=NULL, ...)
}

\arguments{
  \item{inSet}{The input data set as an \code{\link{AffymetrixCelSet}}.}
  \item{outSet}{The output data set as an \code{\link{AffymetrixCelSet}}.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCdf} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFullName.TransformReport]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getInputDataSet.TransformReport]{getInputDataSet}} \tab Gets the source data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.TransformReport]{getName}} \tab Gets the name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getOutputDataSet.TransformReport]{getOutputDataSet}} \tab Gets the transformed data set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.TransformReport]{getPath}} \tab Gets the path of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getTags.TransformReport]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{getUnitNamesFile} \tab  -\cr
 \tab \code{getUnitTypesFile} \tab  -\cr
 \tab \code{getYY} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{plotXYCurve} \tab  -\cr
 \tab \code{plotXYCurveLog2} \tab  -\cr
 \tab \code{seq} \tab  -\cr
 \tab \code{writeImageCombined} \tab  -\cr
 \tab \code{writeImages} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
