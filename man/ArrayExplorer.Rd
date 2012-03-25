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

 \usage{ArrayExplorer(csTuple=NULL, version=c("3.4", "3.0"), ...)}

 \arguments{
   \item{csTuple}{An \code{\link{AffymetrixCelSet}} object.}
   \item{version}{The version of the Explorer HTML/Javascript generated/used.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addColorMap} \tab  -\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getColorMaps} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getDataSet.ArrayExplorer]{getDataSet}} \tab Gets the data set.\cr
 \tab \code{getSetTuple} \tab  -\cr
 \tab \code{nbrOfChipTypes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.ArrayExplorer]{process}} \tab Generates image files, scripts and dynamic pages for the explorer.\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:setArrays.ArrayExplorer]{setArrays}} \tab Sets the arrays.\cr
 \tab \code{setColorMaps} \tab  -\cr
}


 \bold{Methods inherited from Explorer}:\cr
addIncludes, addIndexFile, display, getAlias, getArraysOfInput, getAsteriskTags, getFullName, getIncludePath, getMainPath, getName, getNameOfInput, getNames, getParallelSafe, getPath, getReportPathPattern, getSampleLayerPrefix, getSubname, getTags, getTagsOfInput, getTemplatePath, getVersion, nbrOfArrays, process, setAlias, setArrays, setParallelSafe, setReportPathPattern, setSubname, setup, splitByReportPathPattern, updateSetupExplorerFile

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
