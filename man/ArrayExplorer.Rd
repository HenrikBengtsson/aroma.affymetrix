%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ArrayExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ArrayExplorer}
\docType{class}
\alias{ArrayExplorer}


\title{The ArrayExplorer class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class ArrayExplorer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{Explorer}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ArrayExplorer}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{ArrayExplorer}\cr
extends \emph{\link[aroma.core]{Explorer}}\cr


}

\usage{
ArrayExplorer(csTuple=NULL, ...)
}

\arguments{
  \item{csTuple}{An \code{\link{AffymetrixCelSet}} object.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addColorMap} \tab  -\cr
 \tab \code{getColorMaps} \tab  -\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{getSetTuple} \tab  -\cr
 \tab \code{nbrOfChipTypes} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setArrays} \tab  -\cr
 \tab \code{setColorMaps} \tab  -\cr
}


 \bold{Methods inherited from Explorer}:\cr
addIncludes, addIndexFile, as.character, display, getAlias, getArraysOfInput, getAsteriskTags, getFullName, getIncludePath, getMainPath, getName, getNameOfInput, getNames, getPath, getReportPathPattern, getRootPath, getSampleLayerPrefix, getSubname, getTags, getTagsOfInput, getTemplatePath, getVersion, nbrOfArrays, process, setAlias, setArrays, setReportPathPattern, setSubname, setup, splitByReportPathPattern, updateSetupExplorerFile

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
