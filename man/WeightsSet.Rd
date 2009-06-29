%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  WeightsSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{WeightsSet}
\docType{class}
\alias{WeightsSet}


 \title{The WeightsSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class WeightsSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{WeightsSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{WeightsSet}\cr
extends \link[aroma.affymetrix]{AffymetrixCelSet}\cr



  This class represents probe-level weights.
 }

 \usage{WeightsSet(..., probeModel=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixCelSet}}.}
   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, as, as.AffymetrixCelSet, averageQuantile, byName, convertToUnique, extractAffyBatch, extractMatrix, findByName, fromFiles, fromName, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getChipType, getData, getFullName, getIntensities, getPlatform, getTimestamps, getUnitGroupCellMap, getUnitIntensities, getUnitNamesFile, isDuplicated, justSNPRMA, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, setCdf, update2, writeSgr

 \bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, fromFiles

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitNamesFile

 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
getChipType, getPlatform, nbrOfArrays, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getFullNameTranslator, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setFullName, setFullNamesTranslator, setFullNameTranslator, setName, setTags, update2, validate

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getWeightsSet()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
